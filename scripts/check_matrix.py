"""Read DICOM matrix size for all PCCT cases from the actual DICOM files."""
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ssm_helper import ssm_run

cmd = """python3 << 'PYEOF'
import json, glob, os, struct

base_ws = "/inst/zenith/AppData/Working Storage"
analysts = ["wi_layne.cassidy", "wi_mackenzie.kinney"]

for analyst in analysts:
    for d in sorted(glob.glob(os.path.join(base_ws, analyst, "PT-*"))):
        pt = os.path.basename(d)
        wfiles = glob.glob(os.path.join(d, "**/workitem.json"), recursive=True)
        if not wfiles:
            continue
        try:
            with open(wfiles[0]) as f:
                data = json.load(f)
            iss = data.get("imageSeriesSet", [])
            if not iss:
                continue
            model = iss[0].get("model", "")
            if "NAEOTOM" not in model:
                continue
            series_dir = iss[0].get("seriesLocalFolderName", "")
            dcm_files = [f for f in os.listdir(series_dir) if f.endswith(".dcm")][:1]
            if not dcm_files:
                print(f"{pt}|NO_DCM")
                continue
            with open(os.path.join(series_dir, dcm_files[0]), "rb") as fh:
                raw = fh.read(20000)
            # Rows (0028,0010) and Columns (0028,0011) - implicit VR little endian
            rows = cols = 0
            ri = raw.find(b"\\x28\\x00\\x10\\x00")
            ci = raw.find(b"\\x28\\x00\\x11\\x00")
            if ri >= 0:
                rows = struct.unpack("<H", raw[ri+8:ri+10])[0]
            if ci >= 0:
                cols = struct.unpack("<H", raw[ci+8:ci+10])[0]
            print(f"{pt}|{rows}x{cols}")
        except Exception as e:
            print(f"{pt}|ERROR:{e}")
PYEOF"""

out = ssm_run(cmd, timeout=60)
if out:
    print(f"{'Patient':<25s} {'Matrix'}")
    print("-" * 40)
    for line in out.strip().split("\n"):
        parts = line.split("|")
        if len(parts) >= 2:
            marker = " <-- 1024x1024" if "1024x1024" in parts[1] else ""
            print(f"{parts[0]:<25s} {parts[1]}{marker}")
