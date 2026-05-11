"""Deploy val_helpers to ip3 via chunked base64 over SSM."""
import sys
import os
import json
import subprocess
import time
import glob

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ssm_helper import ssm_run, REGION, PROFILE, get_instance_id

INSTANCE = get_instance_id()

DEST = "/inst/zenith/AppData/erica.freund"


def _aws_args():
    args = ["aws"]
    if PROFILE:
        args += ["--profile", PROFILE]
    args += ["--region", REGION]
    return args

# First create the destination and clear any previous transfer
print("Preparing destination on ip3...")
out = ssm_run(f"mkdir -p '{DEST}' && rm -f /tmp/val_helpers.b64 && echo OK")
print(f"  {out.strip() if out else 'ERROR'}")

# Send each chunk
chunks = sorted(glob.glob("/tmp/vh_chunk_*"))
print(f"\nSending {len(chunks)} chunks...")
for i, chunk_path in enumerate(chunks):
    with open(chunk_path, "r") as f:
        data = f.read()
    # Use append (>>) for all chunks
    cmd = f"cat >> /tmp/val_helpers.b64 << 'CHUNKEOF'\n{data}\nCHUNKEOF"
    params = json.dumps({"commands": [cmd]})
    with open("/tmp/ssm_params.json", "w") as f:
        f.write(params)

    result = subprocess.run(
        _aws_args() + ["ssm", "send-command",
         "--instance-ids", INSTANCE,
         "--document-name", "AWS-RunShellScript",
         "--parameters", "file:///tmp/ssm_params.json",
         "--query", "Command.CommandId", "--output", "text"],
        capture_output=True, text=True
    )
    cmd_id = result.stdout.strip()
    if len(cmd_id) < 36:
        print(f"  Chunk {i+1} FAILED to send: {result.stderr}")
        sys.exit(1)

    # Wait for completion
    time.sleep(2)
    r2 = subprocess.run(
        _aws_args() + ["ssm", "get-command-invocation",
         "--command-id", cmd_id, "--instance-id", INSTANCE,
         "--output", "json"],
        capture_output=True, text=True
    )
    inv = json.loads(r2.stdout)
    status = inv.get("Status", "Unknown")
    if status != "Success":
        print(f"  Chunk {i+1} FAILED: {status} — {inv.get('StandardErrorContent', '')[:200]}")
        sys.exit(1)
    print(f"  Chunk {i+1}/{len(chunks)} sent")

# Decode and extract
print("\nDecoding and extracting on ip3...")
out = ssm_run(
    f"cd '{DEST}' && "
    f"base64 -d /tmp/val_helpers.b64 > /tmp/val_helpers.tar.gz && "
    f"tar xzf /tmp/val_helpers.tar.gz && "
    f"rm -f /tmp/val_helpers.b64 /tmp/val_helpers.tar.gz && "
    f"echo '=== Deployed ===' && "
    f"ls -la '{DEST}/val_helpers/' && "
    f"echo '=== Files ===' && "
    f"find '{DEST}/val_helpers/' -type f"
)
if out:
    print(out)
else:
    print("ERROR: no output from extract step")
