"""Inspect readings.json structure to find jurisdiction field."""
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ssm_helper import ssm_run

BASE = "/inst/zenith/AppData/Working Storage/wi_layne.cassidy"

cmd = """python3 << 'PYEOF'
import json, glob, os

base = '""" + BASE + """'
# Look at one file to understand structure
f = glob.glob(os.path.join(base, 'PT-119', '**', 'readings.json'), recursive=True)[0]
with open(f) as fh:
    d = json.load(fh)

print("=== TOP KEYS ===")
for k, v in d.items():
    print(f"  {k}: {type(v).__name__}")

print()
print("=== _metadata ===")
meta = d.get('_metadata', {})
print(json.dumps(meta, indent=2, default=str)[:3000])

print()
print("=== tags ===")
print(json.dumps(d.get('tags', {}), indent=2, default=str)[:1000])

print()
print("=== name ===")
print(d.get('name', ''))
PYEOF"""

out = ssm_run(cmd)
if out:
    print(out)
