"""Run SNR on all PCCT and EID cases from layne/mackenzie on ip3.

Run from PCCT/scripts/ directory: `python run_snr_all.py`
"""
import sys
import os
import json
import subprocess
import time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ssm_helper import ssm_run

DEST = "/inst/zenith/AppData/erica.freund"
VH = f"{DEST}/val_helpers"
PY = f"{DEST}/venv/bin/python3"

# All workitems from the workflow spreadsheet (layne/mackenzie), current as of 2026-05-07.
# PT-152_NAEOTOM Alpha (wi-bb33e995) and PT-163 are tracked elsewhere -- excluded here.
ALL_WIS = [
    # PCCT
    ("PT-116_Bv_Bv44u_4", "wi-851bbb0b", "wi_layne.cassidy", "PCCT"),
    ("PT-119", "wi-173e25a1", "wi_layne.cassidy", "PCCT"),
    ("PT-120", "wi-9dd5c6c0", "wi_mackenzie.kinney", "PCCT"),
    ("PT-121", "wi-55f45c1e", "wi_layne.cassidy", "PCCT"),
    ("PT-124", "wi-d28f293f", "wi_mackenzie.kinney", "PCCT"),
    ("PT-125", "wi-d3cf471f", "wi_mackenzie.kinney", "PCCT"),
    ("PT-127", "wi-5015b5e1", "wi_layne.cassidy", "PCCT"),
    ("PT-129", "wi-aff6fdcf", "wi_layne.cassidy", "PCCT"),
    ("PT-130", "wi-6f5671bb", "wi_layne.cassidy", "PCCT"),
    ("PT-131", "wi-f5003da9", "wi_layne.cassidy", "PCCT"),
    ("PT-133", "wi-b3859f2c", "wi_layne.cassidy", "PCCT"),
    ("PT-135", "wi-3110d645", "wi_layne.cassidy", "PCCT"),
    ("PT-136", "wi-36b184bc", "wi_mackenzie.kinney", "PCCT"),
    ("PT-138", "wi-b6c7c5d4", "wi_layne.cassidy", "PCCT"),
    ("PT-140", "wi-313231c4", "wi_mackenzie.kinney", "PCCT"),
    ("PT-142", "wi-b0dd5ac3", "wi_layne.cassidy", "PCCT"),
    ("PT-143", "wi-d3a215cd", "wi_layne.cassidy", "PCCT"),
    ("PT-149", "wi-75c06b7e", "wi_layne.cassidy", "PCCT"),
    ("PT-150", "wi-bbb9f25c", "wi_mackenzie.kinney", "PCCT"),
    ("PT-153", "wi-48024a31", "wi_layne.cassidy", "PCCT"),
    ("PT-155", "wi-dca4f60a", "wi_layne.cassidy", "PCCT"),
    ("PT-156", "wi-0d793397", "wi_mackenzie.kinney", "PCCT"),
    ("PT-157", "wi-2e15cb65", "wi_mackenzie.kinney", "PCCT"),
    ("PT-158", "wi-8304e48f", "wi_mackenzie.kinney", "PCCT"),
    ("PT-159", "wi-647edebc", "wi_layne.cassidy", "PCCT"),
    ("PT-161", "wi-76afa109", "wi_layne.cassidy", "PCCT"),
    ("PT-162", "wi-a18e64c6", "wi_layne.cassidy", "PCCT"),
    ("PT-165", "wi-7a2415a8", "wi_layne.cassidy", "PCCT"),
    # EID
    ("PT-116", "wi-ea2750f8", "wi_layne.cassidy", "EID"),
    ("PT-119", "wi-06c5f433", "wi_layne.cassidy", "EID"),
    ("PT-124", "wi-769dc6f1", "wi_mackenzie.kinney", "EID"),
    ("PT-125", "wi-36f6fb1f", "wi_mackenzie.kinney", "EID"),
    ("PT-127", "wi-60ad8e68", "wi_layne.cassidy", "EID"),
    ("PT-129", "wi-8622f0aa", "wi_layne.cassidy", "EID"),
    ("PT-130", "wi-d80e8e1e", "wi_layne.cassidy", "EID"),
    ("PT-135", "wi-58a50c3f", "wi_layne.cassidy", "EID"),
    ("PT-136", "wi-45141053", "wi_mackenzie.kinney", "EID"),
    ("PT-138", "wi-4473f084", "wi_layne.cassidy", "EID"),
    ("PT-140", "wi-9c1158eb", "wi_mackenzie.kinney", "EID"),
    ("PT-142", "wi-e5b7995e", "wi_layne.cassidy", "EID"),
    ("PT-150", "wi-378daa95", "wi_mackenzie.kinney", "EID"),
    ("PT-153", "wi-78ed1e81", "wi_layne.cassidy", "EID"),
    ("PT-155", "wi-506dd6bf", "wi_layne.cassidy", "EID"),
    ("PT-156", "wi-4024e298", "wi_mackenzie.kinney", "EID"),
    ("PT-157", "wi-ea76223b", "wi_mackenzie.kinney", "EID"),
    ("PT-158", "wi-21a4ec32", "wi_mackenzie.kinney", "EID"),
    ("PT-161", "wi-3904013a", "wi_mackenzie.kinney", "EID"),
    ("PT-162", "wi-1b014db2", "wi_layne.cassidy", "EID"),
    ("PT-165", "wi-27008d2e", "wi_layne.cassidy", "EID"),
]

BASE = "/inst/zenith/AppData/Working Storage"

# Step 1: Generate readings list with scan type tag
print("Step 1: Generating readings list on ip3...")
entries = []
for pid, wi, analyst, scan_type in ALL_WIS:
    entries.append(f"{BASE}/{analyst}/{pid}/{wi}")

# Build remote script to find readings.json for each workitem
lines = []
for pid, wi, analyst, scan_type in ALL_WIS:
    wi_dir = f"{BASE}/{analyst}/{pid}/{wi}"
    lines.append(f'rj=$(find "{wi_dir}" -name "readings.json" -print -quit 2>/dev/null); '
                 f'if [ -n "$rj" ]; then echo "{scan_type}|{pid}|$rj"; '
                 f'else echo "{scan_type}|{pid}|NONE"; fi')

cmd_gen = " && ".join(lines)
out = ssm_run(cmd_gen, timeout=60)
if not out:
    print("No output from readings list generation")
    sys.exit(1)

# Write separate readings lists for PCCT and EID
pcct_readings = []
eid_readings = []
for line in out.strip().split("\n"):
    parts = line.split("|")
    if len(parts) < 3 or parts[2] == "NONE":
        print(f"  No readings.json: {parts[0]} {parts[1]}")
        continue
    if parts[0] == "PCCT":
        pcct_readings.append(parts[2])
    else:
        eid_readings.append(parts[2])

# Write to ip3
pcct_list = "\n".join(pcct_readings)
eid_list = "\n".join(eid_readings)
ssm_run(f"echo '{pcct_list}' > '{DEST}/pcct_readings_list.csv'", timeout=10)
ssm_run(f"echo '{eid_list}' > '{DEST}/eid_readings_list.csv'", timeout=10)
print(f"  PCCT: {len(pcct_readings)} readings paths")
print(f"  EID: {len(eid_readings)} readings paths")

# Step 2: Run SNR for PCCT cases
print("\nStep 2: Running SNR for PCCT cases...")
cmd_pcct = f"""
cd '{VH}' && \
mkdir -p '{DEST}/snr_results/logs' && \
PYTHONPATH='{VH}' {PY} compute_snr.py \
    --readings_list_csv '{DEST}/pcct_readings_list.csv' \
    --csv '{DEST}/snr_results/pcct_snr.csv' \
    --edge 10.0 \
    --log-level INFO \
    --log-file '{DEST}/snr_results/logs/snr_pcct.log' \
    2>&1
"""

# Use long polling for NRRD processing
params = json.dumps({"commands": [cmd_pcct]})
with open("/tmp/ssm_params.json", "w") as f:
    f.write(params)
r = subprocess.run(
    ["aws", "ssm", "send-command", "--profile", PROFILE, "--region", REGION,
     "--instance-ids", INSTANCE, "--document-name", "AWS-RunShellScript",
     "--parameters", "file:///tmp/ssm_params.json",
     "--timeout-seconds", "3600",
     "--query", "Command.CommandId", "--output", "text"],
    capture_output=True, text=True
)
cmd_id = r.stdout.strip()
print(f"  Command: {cmd_id}")
for attempt in range(60):
    time.sleep(15)
    r2 = subprocess.run(
        ["aws", "ssm", "get-command-invocation", "--profile", PROFILE, "--region", REGION,
         "--command-id", cmd_id, "--instance-id", INSTANCE, "--output", "json"],
        capture_output=True, text=True
    )
    inv = json.loads(r2.stdout)
    status = inv.get("Status", "Unknown")
    elapsed = (attempt + 1) * 15
    if status in ("Success", "Failed"):
        print(f"  [{elapsed}s] {status}")
        stdout = inv.get("StandardOutputContent", "")
        if stdout:
            # Print just the summary line
            for l in stdout.strip().split("\n"):
                if "Wrote" in l or "Summary" in l or "Success" in l or "Fail" in l:
                    print(f"  {l}")
        stderr = inv.get("StandardErrorContent", "")
        if stderr:
            print(f"  stderr: {stderr[:300]}")
        break
    if attempt % 4 == 0:
        print(f"  [{elapsed}s] processing...")

# Step 3: Run SNR for EID cases
print("\nStep 3: Running SNR for EID cases...")
cmd_eid = f"""
cd '{VH}' && \
PYTHONPATH='{VH}' {PY} compute_snr.py \
    --readings_list_csv '{DEST}/eid_readings_list.csv' \
    --csv '{DEST}/snr_results/eid_snr.csv' \
    --edge 10.0 \
    --log-level INFO \
    --log-file '{DEST}/snr_results/logs/snr_eid.log' \
    2>&1
"""
params = json.dumps({"commands": [cmd_eid]})
with open("/tmp/ssm_params.json", "w") as f:
    f.write(params)
r = subprocess.run(
    ["aws", "ssm", "send-command", "--profile", PROFILE, "--region", REGION,
     "--instance-ids", INSTANCE, "--document-name", "AWS-RunShellScript",
     "--parameters", "file:///tmp/ssm_params.json",
     "--timeout-seconds", "3600",
     "--query", "Command.CommandId", "--output", "text"],
    capture_output=True, text=True
)
cmd_id = r.stdout.strip()
print(f"  Command: {cmd_id}")
for attempt in range(60):
    time.sleep(15)
    r2 = subprocess.run(
        ["aws", "ssm", "get-command-invocation", "--profile", PROFILE, "--region", REGION,
         "--command-id", cmd_id, "--instance-id", INSTANCE, "--output", "json"],
        capture_output=True, text=True
    )
    inv = json.loads(r2.stdout)
    status = inv.get("Status", "Unknown")
    elapsed = (attempt + 1) * 15
    if status in ("Success", "Failed"):
        print(f"  [{elapsed}s] {status}")
        stdout = inv.get("StandardOutputContent", "")
        if stdout:
            for l in stdout.strip().split("\n"):
                if "Wrote" in l or "Summary" in l or "Success" in l or "Fail" in l:
                    print(f"  {l}")
        stderr = inv.get("StandardErrorContent", "")
        if stderr:
            print(f"  stderr: {stderr[:300]}")
        break
    if attempt % 4 == 0:
        print(f"  [{elapsed}s] processing...")

# Step 4: Retrieve results
import base64

PCCT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))
LOCAL_RESULTS = os.path.join(PCCT_ROOT, "gate_results")
os.makedirs(LOCAL_RESULTS, exist_ok=True)

print("\nStep 4: Retrieving results...")
for remote_name, local_name in [("pcct_snr.csv", "snr_pcct.csv"), ("eid_snr.csv", "snr_eid.csv")]:
    out = ssm_run(f"cat '{DEST}/snr_results/{remote_name}' 2>/dev/null", timeout=15)
    if out:
        with open(os.path.join(LOCAL_RESULTS, local_name), "w") as f:
            f.write(out)
        lines = out.strip().split("\n")
        print(f"  {local_name}: {len(lines)-1} cases")
    else:
        print(f"  {local_name}: not available")

print("\nDone.")
