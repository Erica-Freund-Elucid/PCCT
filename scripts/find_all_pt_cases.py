"""Find ALL PT cases across all analyst folders on ip3."""
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ssm_helper import ssm_run

cmd = """python3 << 'PYEOF'
import os, glob

base = "/inst/zenith/AppData/Working Storage"
# List all analyst directories
analysts = sorted([d for d in os.listdir(base) if os.path.isdir(os.path.join(base, d)) and d.startswith("wi_")])
print("ANALYSTS:" + "|".join(analysts))

for analyst in analysts:
    apath = os.path.join(base, analyst)
    pt_dirs = sorted(glob.glob(os.path.join(apath, "PT-*")))
    if pt_dirs:
        for d in pt_dirs:
            print(f"PT|{analyst}|{os.path.basename(d)}")
PYEOF"""

out = ssm_run(cmd, timeout=30)
if not out:
    print("No output")
    sys.exit(1)

lines = out.strip().split("\n")
analysts_line = [l for l in lines if l.startswith("ANALYSTS:")]
if analysts_line:
    analysts = analysts_line[0].split(":")[1].split("|")
    print(f"Analyst directories searched: {len(analysts)}")
    for a in analysts:
        print(f"  {a}")
    print()

pt_lines = [l for l in lines if l.startswith("PT|")]
# Group by analyst
from collections import defaultdict
by_analyst = defaultdict(list)
all_pts = set()
for l in pt_lines:
    _, analyst, pt = l.split("|")
    by_analyst[analyst].append(pt)
    all_pts.add(pt.split("_Bv_")[0] if "_Bv_" in pt else pt)

known = {"wi_layne.cassidy", "wi_mackenzie.kinney"}
new_analysts = set(by_analyst.keys()) - known

print(f"Total PT cases found: {len(pt_lines)}")
print(f"Unique patients: {len(all_pts)}")
print()

for analyst in sorted(by_analyst.keys()):
    pts = by_analyst[analyst]
    tag = " *** NEW ***" if analyst in new_analysts else ""
    print(f"{analyst}{tag}: {len(pts)} PT cases")
    for pt in pts:
        print(f"  {pt}")
    print()
