"""Check which PT cases have Derived folders with Volume.nrrd + Aorta.nrrd on ip3."""
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ssm_helper import ssm_run

cmd = """python3 << 'PYEOF'
import json, glob, os
from pathlib import Path

base_ws = "/inst/zenith/AppData/Working Storage"
analysts = ["wi_layne.cassidy", "wi_mackenzie.kinney"]
readings_paths = []

for analyst in analysts:
    for d in sorted(glob.glob(os.path.join(base_ws, analyst, "PT-*"))):
        pt = os.path.basename(d)
        wfiles = glob.glob(os.path.join(d, "**/workitem.json"), recursive=True)
        if not wfiles:
            print(f"{analyst}|{pt}|NO_WORKITEM||")
            continue
        try:
            with open(wfiles[0]) as f:
                data = json.load(f)
            iss = data.get("imageSeriesSet", [])
            if not iss:
                print(f"{analyst}|{pt}|NO_SERIES||")
                continue
            series_str = iss[0].get("seriesLocalFolderName", "")
            model = iss[0].get("model", "")
            series_path = Path(series_str)
            parts = list(series_path.parts)
            idx = parts.index("Images")
            parts.insert(idx, "Derived")
            derived = Path(*parts)
            vol = derived / "Volume.nrrd"
            aorta = derived / "Aorta.nrrd"
            has_vol = vol.exists()
            has_aorta = aorta.exists()
            if has_vol and has_aorta:
                status = "READY"
            elif has_vol:
                status = "NO_AORTA"
            elif derived.exists():
                status = "NO_VOLUME"
            else:
                status = "NO_DERIVED"
            # Also find readings.json for the readings list
            rfiles = glob.glob(os.path.join(d, "**/readings.json"), recursive=True)
            rpath = rfiles[0] if rfiles else ""
            print(f"{analyst}|{pt}|{model}|{status}|{derived}|{rpath}")
        except Exception as e:
            print(f"{analyst}|{pt}||ERROR:{e}||")
PYEOF"""

out = ssm_run(cmd)
if out:
    print(f"{'Analyst':<20s} {'Patient':<22s} {'Model':<18s} {'Status':<12s}")
    print("-" * 80)
    ready = 0
    total = 0
    for line in out.strip().split("\n"):
        parts = line.split("|")
        if len(parts) >= 4:
            total += 1
            if parts[3] == "READY":
                ready += 1
            print(f"{parts[0]:<20s} {parts[1]:<22s} {parts[2]:<18s} {parts[3]:<12s}")
    print(f"\nReady for SNR: {ready}/{total}")
