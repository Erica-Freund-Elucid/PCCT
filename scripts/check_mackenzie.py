"""Check all mackenzie.kinney workitems with scanner info and summary status."""
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ssm_helper import ssm_run

cmd = """python3 << 'PYEOF'
import json, glob, os

base = "/inst/zenith/AppData/Working Storage/wi_mackenzie.kinney"
for d in sorted(glob.glob(os.path.join(base, "PT-*"))):
    pt = os.path.basename(d)
    # Check each workitem subdir
    wi_dirs = sorted(glob.glob(os.path.join(d, "wi-*")))
    if not wi_dirs:
        wi_dirs = [d]
    for wi_dir in wi_dirs:
        wi_id = os.path.basename(wi_dir)
        wfile = os.path.join(wi_dir, "workitem.json")
        if not os.path.exists(wfile):
            continue
        try:
            with open(wfile) as f:
                data = json.load(f)
            iss = data.get("imageSeriesSet", [])
            model = iss[0].get("model", "") if iss else ""
            kernel = iss[0].get("convolutionKernel", "") if iss else ""
            cat = "EID" if "Force" in model else ("PCCT" if "NAEOTOM" in model else "?")
            csvs = glob.glob(os.path.join(wi_dir, "workitem_summary*.csv"))
            latest = ""
            if csvs:
                latest = os.path.basename(sorted(csvs, key=os.path.getmtime)[-1])
            print(f"{pt}|{wi_id}|{cat}|{model}|{kernel}|{len(csvs)}|{latest}")
        except Exception as e:
            print(f"{pt}|{wi_id}|ERROR|{e}|||")
PYEOF"""

out = ssm_run(cmd, timeout=30)
if out:
    print(f"{'Patient':<15s} {'Workitem':<15s} {'Cat':<5s} {'Model':<18s} {'Kernel':<10s} {'CSVs':<5s} {'Latest'}")
    print("-" * 100)
    for line in out.strip().split("\n"):
        parts = line.split("|")
        if len(parts) >= 7:
            print(f"{parts[0]:<15s} {parts[1]:<15s} {parts[2]:<5s} {parts[3]:<18s} {parts[4]:<10s} {parts[5]:<5s} {parts[6]}")
