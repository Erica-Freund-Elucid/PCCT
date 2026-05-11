"""Find all PT cases across analyst working directories on ip3."""
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ssm_helper import ssm_run

WS_BASE = "/inst/zenith/AppData/Working Storage"

# List all PT directories in both analyst folders, with scanner info
cmd = """python3 << 'PYEOF'
import json, glob, os

base = '/inst/zenith/AppData/Working Storage'
analysts = ['wi_layne.cassidy', 'wi_mackenzie.kinney']

for analyst in analysts:
    apath = os.path.join(base, analyst)
    if not os.path.isdir(apath):
        print(f"NOT_FOUND|{analyst}")
        continue
    for d in sorted(glob.glob(os.path.join(apath, 'PT-*'))):
        pt = os.path.basename(d)
        # Get scanner info from workitem.json
        wfiles = glob.glob(os.path.join(d, '**', 'workitem.json'), recursive=True)
        make = model = kernel = kvp = thickness = ''
        has_summary = False
        if wfiles:
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
            except:
                pass
        # Check for summary CSVs
        summaries = glob.glob(os.path.join(d, '**', 'workitem_summary*.csv'), recursive=True)
        has_summary = len(summaries) > 0
        cat = 'EID' if 'Force' in model or 'force' in model else ('PCCT' if 'NAEOTOM' in model or 'naeotom' in model.lower() else 'UNKNOWN')
        print(f"{analyst}|{pt}|{cat}|{make}|{model}|{kernel}|{kvp}|{thickness}|{has_summary}|{len(summaries)}")
PYEOF"""

out = ssm_run(cmd)
if out:
    print(f"{'Analyst':<22s} {'Patient':<25s} {'Cat':<6s} {'Model':<20s} {'Kernel':<12s} {'kVp':<5s} {'Thick':<8s} {'Summary CSVs'}")
    print("-" * 120)
    eid_count = 0
    pcct_count = 0
    for line in out.strip().split("\n"):
        parts = line.split("|")
        if len(parts) >= 10:
            analyst, pt, cat, make, model, kernel, kvp, thick, has_sum, n_sum = parts
            print(f"{analyst:<22s} {pt:<25s} {cat:<6s} {model:<20s} {kernel:<12s} {kvp:<5s} {thick:<8s} {n_sum}")
            if cat == 'EID':
                eid_count += 1
            elif cat == 'PCCT':
                pcct_count += 1
        else:
            print(line)
    print(f"\nTotal: {eid_count} EID, {pcct_count} PCCT")
