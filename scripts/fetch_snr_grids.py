"""Fetch SNR verification grid PNGs from ip3 via chunked tar transfer."""
import sys
import os
import json
import subprocess
import time
import base64
import tarfile
import io

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ssm_helper import REGION, PROFILE, get_instance_id

INSTANCE = get_instance_id()

PCCT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))
REMOTE_DIR = "/inst/zenith/AppData/erica.freund/snr_results/SNR_verification_grids"
LOCAL_DIR = os.path.join(PCCT_ROOT, "gate_results", "SNR_verification_grids")

os.makedirs(LOCAL_DIR, exist_ok=True)


def _aws_args():
    args = ["aws"]
    if PROFILE:
        args += ["--profile", PROFILE]
    args += ["--region", REGION]
    return args


def ssm_send_and_wait(cmd, timeout_polls=60, poll_interval=10):
    """Send SSM command and poll until complete."""
    import tempfile
    params = json.dumps({"commands": [cmd]})
    with tempfile.NamedTemporaryFile("w", suffix=".json", delete=False) as f:
        f.write(params)
        params_path = f.name
    try:
        r = subprocess.run(
            _aws_args() + ["ssm", "send-command",
             "--instance-ids", INSTANCE, "--document-name", "AWS-RunShellScript",
             "--parameters", f"file://{params_path}",
             "--timeout-seconds", "3600",
             "--query", "Command.CommandId", "--output", "text"],
            capture_output=True, text=True
        )
    finally:
        try:
            os.unlink(params_path)
        except OSError:
            pass
    cmd_id = r.stdout.strip()
    if len(cmd_id) < 36:
        print(f"  ERROR sending: {r.stderr}")
        return None

    for i in range(timeout_polls):
        time.sleep(poll_interval)
        r2 = subprocess.run(
            _aws_args() + ["ssm", "get-command-invocation",
             "--command-id", cmd_id, "--instance-id", INSTANCE, "--output", "json"],
            capture_output=True, text=True
        )
        inv = json.loads(r2.stdout)
        status = inv.get("Status", "Unknown")
        if status in ("Success", "Failed", "Cancelled", "TimedOut"):
            return inv
    return None


# Step 1: Tar the grids on ip3 and split into base64 chunks
print("Tarring grids on ip3...")
inv = ssm_send_and_wait(
    f"cd /inst/zenith/AppData/erica.freund/snr_results && "
    f"tar czf /tmp/snr_grids.tar.gz SNR_verification_grids/ && "
    f"base64 /tmp/snr_grids.tar.gz > /tmp/snr_grids.b64 && "
    f"wc -c /tmp/snr_grids.b64 && "
    f"split -b 20000 /tmp/snr_grids.b64 /tmp/snr_chunk_ && "
    f"ls /tmp/snr_chunk_* | wc -l",
    poll_interval=5
)
if not inv or inv.get("Status") != "Success":
    print(f"  Failed: {inv}")
    sys.exit(1)
print(f"  {inv.get('StandardOutputContent', '').strip()}")

# Step 2: Get chunk list
inv2 = ssm_send_and_wait("ls /tmp/snr_chunk_*", poll_interval=3)
chunks = [c.strip() for c in inv2["StandardOutputContent"].strip().split("\n") if c.strip()]
print(f"\n{len(chunks)} chunks to transfer...")

# Step 3: Download each chunk
all_b64 = ""
for i, chunk_path in enumerate(chunks):
    inv_c = ssm_send_and_wait(f"cat '{chunk_path}'", poll_interval=3)
    if inv_c and inv_c.get("Status") == "Success":
        all_b64 += inv_c["StandardOutputContent"]
        print(f"  [{i+1}/{len(chunks)}] received")
    else:
        print(f"  [{i+1}/{len(chunks)}] FAILED")
        sys.exit(1)

# Step 4: Decode and extract
print("\nDecoding and extracting...")
tar_bytes = base64.b64decode(all_b64.strip())
print(f"  Tar size: {len(tar_bytes)//1024}KB")

tar = tarfile.open(fileobj=io.BytesIO(tar_bytes), mode="r:gz")
count = 0
for member in tar.getmembers():
    if member.isfile() and member.name.endswith(".png"):
        fname = os.path.basename(member.name)
        content = tar.extractfile(member).read()
        with open(os.path.join(LOCAL_DIR, fname), "wb") as f:
            f.write(content)
        print(f"  {fname} ({len(content)//1024}KB)")
        count += 1

# Cleanup remote temp files
ssm_send_and_wait("rm -f /tmp/snr_grids.tar.gz /tmp/snr_grids.b64 /tmp/snr_chunk_*", poll_interval=3)

print(f"\nDone. {count} PNGs saved to {LOCAL_DIR}")
