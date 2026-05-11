"""Get jurisdiction from readings.json for each PT patient on ip3."""
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ssm_helper import ssm_run

BASE = "/inst/zenith/AppData/Working Storage/wi_layne.cassidy"

# First find all readings.json files in PT directories
find_cmd = (
    f"for d in '{BASE}'/PT-*/; do "
    "pt=$(basename \"$d\"); "
    "rj=$(find \"$d\" -name 'readings.json' -print -quit 2>/dev/null); "
    "if [ -n \"$rj\" ]; then echo \"$pt|$rj\"; "
    "else echo \"$pt|NONE\"; fi; done"
)

print("Finding readings.json files...")
output = ssm_run(find_cmd)
if not output:
    print("No output")
    sys.exit(1)

lines = [l.strip() for l in output.strip().split("\n") if l.strip()]

# For each that has a readings.json, extract jurisdiction
pts_with_json = [(l.split("|")[0], l.split("|")[1]) for l in lines if "|NONE" not in l]
pts_without = [l.split("|")[0] for l in lines if "|NONE" in l]

if pts_without:
    print(f"\nNo readings.json found for: {', '.join(pts_without)}")

# Build a command that reads jurisdiction from each readings.json
# Process in one SSM call to avoid 16 round-trips
extract_cmd = "python3 -c \"\nimport json, glob, os\nbase = '{}'\nfor d in sorted(glob.glob(os.path.join(base, 'PT-*'))):\n    pt = os.path.basename(d)\n    rfiles = glob.glob(os.path.join(d, '**', 'readings.json'), recursive=True)\n    if not rfiles:\n        print(f'{{pt}}|NO_FILE')\n        continue\n    try:\n        with open(rfiles[0]) as f:\n            data = json.load(f)\n        # Try common key names\n        j = None\n        if isinstance(data, dict):\n            for key in ['jurisdiction', 'Jurisdiction', 'site', 'Site', 'institution', 'Institution']:\n                if key in data:\n                    j = data[key]\n                    break\n            # Maybe nested in metadata or patient info\n            if j is None:\n                for k, v in data.items():\n                    if isinstance(v, dict):\n                        for key in ['jurisdiction', 'Jurisdiction', 'site', 'Site']:\n                            if key in v:\n                                j = v[key]\n                                break\n                    if j is not None:\n                        break\n        if j is None:\n            # Print top-level keys for debugging\n            keys = list(data.keys())[:10] if isinstance(data, dict) else ['NOT_DICT']\n            print(f'{{pt}}|KEYS:{{keys}}')\n        else:\n            print(f'{{pt}}|{{j}}')\n    except Exception as e:\n        print(f'{{pt}}|ERROR:{{e}}')\n\"".format(BASE)

print("\nExtracting jurisdiction info...")
output2 = ssm_run(extract_cmd)
if not output2:
    print("No output from extraction")
    sys.exit(1)

print(f"\n{'Patient':<25s} {'Jurisdiction'}")
print("-" * 60)
for line in sorted(output2.strip().split("\n")):
    if "|" in line:
        pt, val = line.split("|", 1)
        print(f"{pt:<25s} {val}")
