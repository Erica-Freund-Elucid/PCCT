"""Full inventory of all PT cases on ip3: scanner, matrix, derived status, summary CSVs."""
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ssm_helper import ssm_run

cmd = """python3 << 'PYEOF'
import json, glob, os, struct
from pathlib import Path

base = "/inst/zenith/AppData/Working Storage"
analysts = sorted([d for d in os.listdir(base)
                   if os.path.isdir(os.path.join(base, d)) and d.startswith("wi_")])

for analyst in analysts:
    for d in sorted(glob.glob(os.path.join(base, analyst, "PT-*"))):
        pt = os.path.basename(d)
        wfiles = glob.glob(os.path.join(d, "**/workitem.json"), recursive=True)
        if not wfiles:
            print(f"{analyst}|{pt}|||||||NO_WORKITEM")
            continue
        try:
            with open(wfiles[0]) as f:
                data = json.load(f)
            iss = data.get("imageSeriesSet", [])
            if not iss:
                print(f"{analyst}|{pt}|||||||NO_SERIES")
                continue
            s = iss[0]
            model = s.get("model", "")
            kernel = s.get("convolutionKernel", "")
            kvp = s.get("kVp", "")
            thickness = s.get("sliceThickness", "")
            series_dir = s.get("seriesLocalFolderName", "")
            cat = "EID" if "Force" in model else ("PCCT" if "NAEOTOM" in model else "OTHER")

            # Matrix from DICOM
            matrix = ""
            try:
                dcm_files = [f for f in os.listdir(series_dir) if f.endswith(".dcm")][:1]
                if dcm_files:
                    with open(os.path.join(series_dir, dcm_files[0]), "rb") as fh:
                        raw = fh.read(20000)
                    ri = raw.find(b"\x28\x00\x10\x00")
                    ci = raw.find(b"\x28\x00\x11\x00")
                    rows = struct.unpack("<H", raw[ri+8:ri+10])[0] if ri >= 0 else 0
                    cols = struct.unpack("<H", raw[ci+8:ci+10])[0] if ci >= 0 else 0
                    matrix = f"{rows}x{cols}"
            except:
                matrix = "ERR"

            # Derived folder status
            try:
                parts = list(Path(series_dir).parts)
                idx = parts.index("Images")
                parts.insert(idx, "Derived")
                derived = Path(*parts)
                has_vol = (derived / "Volume.nrrd").exists()
                has_aorta = (derived / "Aorta.nrrd").exists()
                if has_vol and has_aorta:
                    derived_status = "READY"
                elif has_vol:
                    derived_status = "NO_AORTA"
                elif derived.exists():
                    derived_status = "NO_VOL"
                else:
                    derived_status = "NO_DERIVED"
            except:
                derived_status = "ERR"

            # Summary CSVs
            summaries = glob.glob(os.path.join(d, "**/workitem_summary*.csv"), recursive=True)
            n_sum = len(summaries)

            print(f"{analyst}|{pt}|{cat}|{model}|{kernel}|{kvp}|{matrix}|{derived_status}|{n_sum}")
        except Exception as e:
            print(f"{analyst}|{pt}|||||||ERROR:{e}")
PYEOF"""

out = ssm_run(cmd, timeout=120)
if not out:
    print("No output")
    sys.exit(1)

lines = out.strip().split("\n")

# Parse and display
print(f"{'Analyst':<28s} {'Patient':<22s} {'Cat':<5s} {'Model':<16s} {'Kernel':<10s} {'kVp':<5s} {'Matrix':<10s} {'Derived':<10s} {'CSVs'}")
print("-" * 140)

from collections import defaultdict
patients = defaultdict(list)  # pt_norm -> list of entries

for line in lines:
    parts = line.split("|")
    if len(parts) < 9:
        continue
    analyst, pt, cat, model, kernel, kvp, matrix, derived, csvs = parts[:9]
    print(f"{analyst:<28s} {pt:<22s} {cat:<5s} {model:<16s} {kernel:<10s} {kvp:<5s} {matrix:<10s} {derived:<10s} {csvs}")
    pt_norm = pt.split("_Bv_")[0] if "_Bv_" in pt else pt
    patients[pt_norm].append({
        "analyst": analyst, "pt": pt, "cat": cat, "model": model,
        "matrix": matrix, "derived": derived, "csvs": int(csvs) if csvs.isdigit() else 0
    })

# Summary
print(f"\n{'='*60}")
print(f"SUMMARY")
print(f"{'='*60}")
all_cats = defaultdict(int)
ready_for_snr = 0
has_summary = 0
unique_1024 = set()
for pt_norm, entries in patients.items():
    for e in entries:
        all_cats[e["cat"]] += 1
        if e["derived"] == "READY":
            ready_for_snr += 1
        if e["csvs"] > 0:
            has_summary += 1
        if "1024" in e["matrix"]:
            unique_1024.add(pt_norm)

print(f"Unique patients: {len(patients)}")
print(f"Total workitems: {len(lines)}")
print(f"PCCT: {all_cats['PCCT']}, EID: {all_cats['EID']}, Other: {all_cats.get('OTHER',0)}")
print(f"Ready for SNR (Derived + Aorta): {ready_for_snr}")
print(f"With summary CSVs: {has_summary}")
if unique_1024:
    print(f"1024x1024 patients: {sorted(unique_1024)}")
else:
    print("1024x1024 patients: none found")
