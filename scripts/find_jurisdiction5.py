"""Extract scanner info from workitem.json for each PT patient to infer site."""
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ssm_helper import ssm_run

BASE = "/inst/zenith/AppData/Working Storage/wi_layne.cassidy"

cmd = """python3 << 'PYEOF'
import json, glob, os

base = '""" + BASE + """'
for d in sorted(glob.glob(os.path.join(base, 'PT-*'))):
    pt = os.path.basename(d)
    wfiles = glob.glob(os.path.join(d, '**', 'workitem.json'), recursive=True)
    if not wfiles:
        print(f"{pt}|NO_WORKITEM")
        continue
    try:
        with open(wfiles[0]) as f:
            data = json.load(f)
        iss = data.get('imageSeriesSet', [])
        if iss:
            s = iss[0]
            make = s.get('make', '')
            model = s.get('model', '')
            kernel = s.get('convolutionKernel', '')
            kvp = s.get('kVp', '')
            thickness = s.get('sliceThickness', '')
            anatomy = s.get('anatomy', '')
            print(f"{pt}|{make}|{model}|{kernel}|{kvp}|{thickness}")
        else:
            print(f"{pt}|NO_SERIES")
    except Exception as e:
        print(f"{pt}|ERROR:{e}")
PYEOF"""

out = ssm_run(cmd)
if out:
    print(f"{'Patient':<25s} {'Make':<25s} {'Model':<20s} {'Kernel':<15s} {'kVp':<6s} {'Thickness'}")
    print("-" * 110)
    for line in sorted(out.strip().split("\n")):
        parts = line.split("|")
        if len(parts) >= 6:
            print(f"{parts[0]:<25s} {parts[1]:<25s} {parts[2]:<20s} {parts[3]:<15s} {parts[4]:<6s} {parts[5]}")
        else:
            print(line)
