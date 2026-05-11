"""Refresh all PCCT and EID summary CSVs based on latest workflow spreadsheet.

Run from PCCT/scripts/ directory: `python refresh_all_csvs.py`
"""
import sys
import os
import base64
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ssm_helper import ssm_run

PCCT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))
BASE = "/inst/zenith/AppData/Working Storage"
LOCAL_PCCT = os.path.join(PCCT_ROOT, "workitem_summaries", "PCCT")
LOCAL_EID = os.path.join(PCCT_ROOT, "workitem_summaries", "EID")

# From the workflow spreadsheet (current as of 2026-05-07).
# PT-152_NAEOTOM Alpha (wi-bb33e995) is in the spreadsheet but not yet on the
# instance; PT-163 (wi-5622d81e) held for separate 1024x1024 evaluation.
PCCT_WIS = [
    ("PT-116_Bv_Bv44u_4", "wi-851bbb0b", "wi_layne.cassidy"),
    ("PT-119", "wi-173e25a1", "wi_layne.cassidy"),
    ("PT-120", "wi-9dd5c6c0", "wi_mackenzie.kinney"),
    ("PT-121", "wi-55f45c1e", "wi_layne.cassidy"),
    ("PT-124", "wi-d28f293f", "wi_mackenzie.kinney"),
    ("PT-125", "wi-d3cf471f", "wi_mackenzie.kinney"),
    ("PT-127", "wi-5015b5e1", "wi_layne.cassidy"),
    ("PT-129", "wi-aff6fdcf", "wi_layne.cassidy"),
    ("PT-130", "wi-6f5671bb", "wi_layne.cassidy"),
    ("PT-131", "wi-f5003da9", "wi_layne.cassidy"),
    ("PT-133", "wi-b3859f2c", "wi_layne.cassidy"),
    ("PT-135", "wi-3110d645", "wi_layne.cassidy"),
    ("PT-136", "wi-36b184bc", "wi_mackenzie.kinney"),
    ("PT-138", "wi-b6c7c5d4", "wi_layne.cassidy"),
    ("PT-140", "wi-313231c4", "wi_mackenzie.kinney"),
    ("PT-142", "wi-b0dd5ac3", "wi_layne.cassidy"),
    ("PT-143", "wi-d3a215cd", "wi_layne.cassidy"),
    ("PT-149", "wi-75c06b7e", "wi_layne.cassidy"),
    ("PT-150", "wi-bbb9f25c", "wi_mackenzie.kinney"),
    ("PT-152_NAEOTOM Alpha", "wi-bb33e995", "wi_mackenzie.kinney"),
    ("PT-153", "wi-48024a31", "wi_layne.cassidy"),
    ("PT-155", "wi-dca4f60a", "wi_layne.cassidy"),
    ("PT-156", "wi-0d793397", "wi_mackenzie.kinney"),
    ("PT-157", "wi-2e15cb65", "wi_mackenzie.kinney"),
    ("PT-158", "wi-8304e48f", "wi_mackenzie.kinney"),
    ("PT-159", "wi-647edebc", "wi_layne.cassidy"),
    ("PT-161", "wi-76afa109", "wi_layne.cassidy"),
    ("PT-162", "wi-a18e64c6", "wi_layne.cassidy"),
    ("PT-165", "wi-7a2415a8", "wi_layne.cassidy"),
]

EID_WIS = [
    ("PT-116", "wi-ea2750f8", "wi_layne.cassidy"),
    ("PT-119", "wi-06c5f433", "wi_layne.cassidy"),
    ("PT-124", "wi-769dc6f1", "wi_mackenzie.kinney"),
    ("PT-125", "wi-36f6fb1f", "wi_mackenzie.kinney"),
    ("PT-127", "wi-60ad8e68", "wi_layne.cassidy"),
    ("PT-129", "wi-8622f0aa", "wi_layne.cassidy"),
    ("PT-130", "wi-d80e8e1e", "wi_layne.cassidy"),
    ("PT-135", "wi-58a50c3f", "wi_layne.cassidy"),
    ("PT-136", "wi-45141053", "wi_mackenzie.kinney"),
    ("PT-138", "wi-4473f084", "wi_layne.cassidy"),
    ("PT-140", "wi-9c1158eb", "wi_mackenzie.kinney"),
    ("PT-142", "wi-e5b7995e", "wi_layne.cassidy"),
    ("PT-150", "wi-378daa95", "wi_mackenzie.kinney"),
    ("PT-153", "wi-78ed1e81", "wi_layne.cassidy"),
    ("PT-155", "wi-506dd6bf", "wi_layne.cassidy"),
    ("PT-156", "wi-4024e298", "wi_mackenzie.kinney"),
    ("PT-157", "wi-ea76223b", "wi_mackenzie.kinney"),
    ("PT-158", "wi-21a4ec32", "wi_mackenzie.kinney"),
    ("PT-161", "wi-3904013a", "wi_mackenzie.kinney"),
    ("PT-162", "wi-1b014db2", "wi_layne.cassidy"),
    ("PT-165", "wi-27008d2e", "wi_layne.cassidy"),
]
# PT-142: use Layne's wi-e5b7995e per analyst preference (summary CSV not yet generated as of 2026-05-07).
# PT-163 EID (wi-526b81de) held for separate 1024x1024 evaluation -- not pulled here.


def fetch_latest_csv(pid, wi, analyst, local_dir):
    """Find latest summary CSV for a workitem and download it."""
    remote_dir = f"{BASE}/{analyst}/{pid}/{wi}"
    cmd = f"ls -t '{remote_dir}'/workitem_summary*.csv 2>/dev/null | head -1"
    out = ssm_run(cmd, timeout=10)
    if not out or not out.strip():
        return None, 0

    remote_path = out.strip()
    remote_size_out = ssm_run(f"wc -c < '{remote_path}'", timeout=10)
    rsize = int(remote_size_out.strip()) if remote_size_out else 0

    # Try single fetch first
    b64 = ssm_run(f"base64 '{remote_path}'", timeout=15)
    data = base64.b64decode(b64.strip()) if b64 else b""

    # Chunked if truncated
    if len(data) < rsize:
        ssm_run(f"base64 '{remote_path}' > /tmp/csv_b64 && split -b 18000 /tmp/csv_b64 /tmp/csv_chunk_", timeout=15)
        chunk_list = ssm_run("ls /tmp/csv_chunk_*", timeout=10)
        chunks = [c.strip() for c in chunk_list.strip().split("\n") if c.strip()]
        all_b64 = ""
        for chunk in chunks:
            part = ssm_run(f"cat '{chunk}'", timeout=10)
            if part:
                all_b64 += part
        data = base64.b64decode(all_b64.strip())
        ssm_run("rm -f /tmp/csv_b64 /tmp/csv_chunk_*", timeout=5)

    # Normalize patient ID for filename
    pid_norm = pid.split("_Bv_")[0] if "_Bv_" in pid else pid
    pid_norm = pid_norm.split("_NAEOTOM")[0] if "_NAEOTOM" in pid_norm else pid_norm
    fname = f"{pid_norm}_{os.path.basename(remote_path)}"
    local_path = os.path.join(local_dir, fname)

    with open(local_path, "wb") as f:
        f.write(data)
    return fname, len(data)


# Clear and rebuild both folders
for d in [LOCAL_PCCT, LOCAL_EID]:
    os.makedirs(d, exist_ok=True)
    for f in os.listdir(d):
        if f.endswith(".csv"):
            os.remove(os.path.join(d, f))

# Fetch PCCT
print(f"=== PCCT ({len(PCCT_WIS)} workitems) ===\n")
pcct_ok = 0
for pid, wi, analyst in PCCT_WIS:
    fname, size = fetch_latest_csv(pid, wi, analyst, LOCAL_PCCT)
    if fname:
        print(f"  {fname} ({size//1024}KB)")
        pcct_ok += 1
    else:
        print(f"  {pid} ({wi}): NO CSV YET")
print(f"\n  {pcct_ok}/{len(PCCT_WIS)} PCCT CSVs fetched")

# Fetch EID
print(f"\n=== EID/CCTA ({len(EID_WIS)} workitems) ===\n")
eid_ok = 0
for pid, wi, analyst in EID_WIS:
    fname, size = fetch_latest_csv(pid, wi, analyst, LOCAL_EID)
    if fname:
        print(f"  {fname} ({size//1024}KB)")
        eid_ok += 1
    else:
        print(f"  {pid} ({wi}): NO CSV YET")
print(f"\n  {eid_ok}/{len(EID_WIS)} EID CSVs fetched")

print(f"\nDone. Ready to run gate analyses.")
