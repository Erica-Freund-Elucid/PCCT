"""Helper to run SSM commands on the ip3-manager instance and retrieve output.

Instance ID is resolved dynamically by Name tag so the script self-heals when
the EC2 instance is recycled. Set IP3_INSTANCE_ID env var to skip the lookup.
Set AWS_PROFILE if you need a non-default profile.
"""
import json
import os
import subprocess
import sys
import tempfile
import time

REGION = os.environ.get("AWS_REGION", "us-east-1")
PROFILE = os.environ.get("AWS_PROFILE")  # None = default profile / SSO
INSTANCE_NAME_FILTER = "*ip3-manager*"

_INSTANCE_CACHE = None


def _aws(cmd_args):
    """Run aws CLI, returning (stdout, stderr, returncode)."""
    cmd = ["aws"]
    if PROFILE:
        cmd += ["--profile", PROFILE]
    cmd += ["--region", REGION] + cmd_args
    r = subprocess.run(cmd, capture_output=True, text=True)
    return r.stdout, r.stderr, r.returncode


def get_instance_id():
    """Resolve the ip3-manager instance ID dynamically (cached for the process)."""
    global _INSTANCE_CACHE
    if _INSTANCE_CACHE:
        return _INSTANCE_CACHE
    if os.environ.get("IP3_INSTANCE_ID"):
        _INSTANCE_CACHE = os.environ["IP3_INSTANCE_ID"]
        return _INSTANCE_CACHE
    out, err, rc = _aws([
        "ssm", "describe-instance-information",
        "--filters", f"Key=tag:Name,Values={INSTANCE_NAME_FILTER}",
        "--query", "InstanceInformationList[?PingStatus=='Online'].InstanceId",
        "--output", "text",
    ])
    if rc != 0 or not out.strip():
        # Fall back to filtering by ComputerName since some instances lack the tag
        out, err, rc = _aws([
            "ssm", "describe-instance-information",
            "--query", "InstanceInformationList[?PingStatus=='Online'].[InstanceId,ComputerName]",
            "--output", "text",
        ])
        if rc != 0:
            print(f"ERROR resolving instance: {err}", file=sys.stderr)
            return None
        for line in out.strip().split("\n"):
            parts = line.split()
            if len(parts) >= 2 and "ip3-manager" in parts[1]:
                _INSTANCE_CACHE = parts[0]
                return _INSTANCE_CACHE
        print("ERROR: no online ip3-manager instance found", file=sys.stderr)
        return None
    _INSTANCE_CACHE = out.strip().split()[0]
    return _INSTANCE_CACHE


def ssm_run(shell_cmd, timeout=30):
    """Run a shell command on ip3-manager via SSM and return stdout text."""
    instance = get_instance_id()
    if not instance:
        return None

    params = json.dumps({"commands": [shell_cmd]})
    with tempfile.NamedTemporaryFile("w", suffix=".json", delete=False) as f:
        f.write(params)
        params_path = f.name
    try:
        out, err, rc = _aws([
            "ssm", "send-command",
            "--instance-ids", instance,
            "--document-name", "AWS-RunShellScript",
            "--parameters", f"file://{params_path}",
            "--query", "Command.CommandId", "--output", "text",
        ])
    finally:
        try:
            os.unlink(params_path)
        except OSError:
            pass
    cmd_id = out.strip()
    if rc != 0 or not cmd_id or len(cmd_id) < 36:
        print(f"ERROR sending command: {err or out}", file=sys.stderr)
        return None

    poll_interval = 5
    max_polls = max(3, int(timeout / poll_interval))
    for _ in range(max_polls):
        time.sleep(poll_interval)
        out, err, rc = _aws([
            "ssm", "get-command-invocation",
            "--command-id", cmd_id,
            "--instance-id", instance,
            "--output", "json",
        ])
        if rc != 0:
            continue
        data = json.loads(out)
        status = data.get("Status", "Unknown")
        if status in ("Success", "Failed", "Cancelled", "TimedOut"):
            if status != "Success":
                print(f"Command status: {status}", file=sys.stderr)
                print(f"Stderr: {data.get('StandardErrorContent', '')}", file=sys.stderr)
            return data.get("StandardOutputContent", "")

    print(f"Command timed out after {timeout}s", file=sys.stderr)
    return None


if __name__ == "__main__":
    BASE = "/inst/zenith/AppData/Working Storage/wi_layne.cassidy"

    cmd = f"""cd '{BASE}' && for d in PT-*/; do
        pt=$(basename "$d")
        latest=$(find "$d" -name 'workitem_summary*.csv' -printf '%T+ %p\\n' 2>/dev/null | sort -r | head -1 | cut -d' ' -f2-)
        if [ -n "$latest" ]; then echo "$pt|$latest"; fi
    done"""

    print("Finding latest summary CSVs per PT patient...")
    output = ssm_run(cmd)
    if not output:
        print("No output from find command")
        sys.exit(1)

    lines = [l.strip() for l in output.strip().split("\n") if l.strip()]
    print(f"Found {len(lines)} PT patients with summary CSVs:\n")
    for line in lines:
        pt, path = line.split("|", 1)
        print(f"  {pt}: {path.split('/')[-1]}")

    import base64

    dest = sys.argv[1] if len(sys.argv) > 1 else "/tmp/pt_summaries"
    os.makedirs(dest, exist_ok=True)

    count = 0
    for line in lines:
        pt, path = line.split("|", 1)
        fname = f"{pt}_{path.split('/')[-1]}"
        b64 = ssm_run(f"base64 '{BASE}/{path}'")
        if not b64:
            print(f"  FAILED: {fname}")
            continue
        content = base64.b64decode(b64.strip())
        out_path = os.path.join(dest, fname)
        with open(out_path, "wb") as f:
            f.write(content)
        print(f"  -> {fname} ({len(content)} bytes)")
        count += 1

    print(f"\n{count} files saved to {dest}")
