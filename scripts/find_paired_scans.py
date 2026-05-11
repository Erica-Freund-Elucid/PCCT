"""Find PT patients with paired PCCT/EID scans across all analyst folders."""
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ssm_helper import ssm_run

cmd = """python3 << 'PYEOF'
import json, glob, os

base = '/inst/zenith/AppData/Working Storage'
analysts = ['wi_layne.cassidy', 'wi_mackenzie.kinney']

# Collect all workitems per PT patient across all analysts
patient_scans = {}  # pt_id -> list of {analyst, wi, model, category, has_summary}

for analyst in analysts:
    apath = os.path.join(base, analyst)
    if not os.path.isdir(apath):
        continue
    for d in sorted(glob.glob(os.path.join(apath, 'PT-*'))):
        pt = os.path.basename(d)
        # Normalize PT ID (strip kernel suffix like _Bv_Bv44u_4)
        pt_norm = pt.split('_Bv_')[0] if '_Bv_' in pt else pt

        # Each PT dir may have multiple workitem subdirs
        wi_dirs = [x for x in glob.glob(os.path.join(d, 'wi-*')) if os.path.isdir(x)]
        if not wi_dirs:
            wi_dirs = [d]  # fallback

        for wi_dir in wi_dirs:
            wi_id = os.path.basename(wi_dir)
            wfile = os.path.join(wi_dir, 'workitem.json')
            if not os.path.exists(wfile):
                # try nested
                wfiles = glob.glob(os.path.join(wi_dir, '**', 'workitem.json'), recursive=True)
                wfile = wfiles[0] if wfiles else None

            model = kernel = ''
            cat = 'UNKNOWN'
            if wfile and os.path.exists(wfile):
                try:
                    with open(wfile) as f:
                        data = json.load(f)
                    iss = data.get('imageSeriesSet', [])
                    if iss:
                        model = iss[0].get('model', '')
                        kernel = iss[0].get('convolutionKernel', '')
                        if 'Force' in model:
                            cat = 'EID'
                        elif 'NAEOTOM' in model:
                            cat = 'PCCT'
                except:
                    pass

            summaries = glob.glob(os.path.join(wi_dir, 'workitem_summary*.csv'))

            if pt_norm not in patient_scans:
                patient_scans[pt_norm] = []
            patient_scans[pt_norm].append({
                'folder': pt,
                'analyst': analyst.replace('wi_', ''),
                'wi': wi_id,
                'model': model,
                'kernel': kernel,
                'cat': cat,
                'has_summary': len(summaries) > 0,
                'n_summaries': len(summaries),
            })

# Report
print("=== ALL PT PATIENTS AND THEIR SCANS ===")
print()
paired = 0
pcct_only = 0
eid_only = 0
for pt in sorted(patient_scans.keys()):
    scans = patient_scans[pt]
    cats = set(s['cat'] for s in scans)
    has_pcct = 'PCCT' in cats
    has_eid = 'EID' in cats
    status = 'PAIRED' if has_pcct and has_eid else ('PCCT_ONLY' if has_pcct else ('EID_ONLY' if has_eid else 'UNKNOWN'))
    if status == 'PAIRED':
        paired += 1
    elif status == 'PCCT_ONLY':
        pcct_only += 1
    elif status == 'EID_ONLY':
        eid_only += 1

    print(f"{pt}|{status}")
    for s in scans:
        sum_status = f"{s['n_summaries']} CSV" if s['has_summary'] else "NO CSV"
        print(f"  {s['analyst']}|{s['folder']}|{s['wi']}|{s['cat']}|{s['model']}|{s['kernel']}|{sum_status}")

print()
print(f"SUMMARY: {paired} paired, {pcct_only} PCCT-only, {eid_only} EID-only, {len(patient_scans)} total patients")
print(f"TARGET: 30 paired patients")
print(f"REMAINING: {max(0, 30 - paired)} pairs needed")
PYEOF"""

out = ssm_run(cmd)
if out:
    print(out)
