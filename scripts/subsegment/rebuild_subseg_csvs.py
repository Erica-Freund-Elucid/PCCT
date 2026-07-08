"""Rebuild sub-segment workitem_summary CSVs from the intersection-pipeline output.

Runs on the compute box (reads /var/tmp/pcct_subseg/out/<PID>/final/{left,right}/).
For each patient it parses common_vessel_sections_volume_summary.txt (per-vessel
Study A = PCCT, Study B = EID volumes) and common_vessel_sections_summary.json
(per-vessel common length), then writes vessel-level CSVs matching the schema the
gate analysis expects, applying the verified mapping:

    LumenVol        = lumen
    LumenAndWallVol = wall            (wallPartition region is LumenAndWall)
    WallVol         = wall - lumen    (wall-only tissue)
    CALCVol/LRNCVol/NonCALCMATXVol = respective plaque regions (direct)
    TotalPlaqueVolume = CALC + NonCALCMATX

Study A -> PCCT CSV, Study B -> EID CSV. Usage:
    python rebuild_subseg_csvs.py <out_root> <pairs.txt> <out_csv_dir>
pairs.txt lines: PID|pcct_dir|eid_dir  (wi = basename of each dir)
"""
import csv, json, os, re, sys

COMPONENTS = {"Wall": "wall", "Lumen": "lumen", "Plaque - CALC": "CALC",
              "Plaque - LRNC": "LRNC", "Plaque - NonCALCMATX": "NonCALCMATX"}
HEADER = ["workitemID", "individualID", "subjectID", "alsoKnownAs", "targetID",
          "bodySite", "applies", "age", "dob", "indication", "sex", "location",
          "level", "Len", "LumenVol", "WallVol", "LumenAndWallVol",
          "CALCVol", "LRNCVol", "NonCALCMATXVol", "TotalPlaqueVolume"]
ROW = re.compile(r"^(\S.*?)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+-?\d+\.\d+\s*$")


def parse_volume_txt(path):
    """Return {component: {vessel: (studyA, studyB)}}."""
    comp = {}
    cur = None
    for line in open(path, encoding="utf-8"):
        s = line.rstrip("\n")
        key = s.strip()
        if key in COMPONENTS:
            cur = COMPONENTS[key]; comp[cur] = {}
            continue
        if cur is None or key.startswith("Total") or line[:4].strip() == "" and line.startswith("    "):
            continue
        m = ROW.match(s)
        if m and not m.group(1).startswith("Total"):
            comp[cur][m.group(1).strip()] = (float(m.group(2)), float(m.group(3)))
    return comp


def vessel_lengths(json_path):
    d = json.load(open(json_path, encoding="utf-8"))
    return ({mv["vessel_name"]: mv.get("total_common_length_mm", "") for mv in d.get("matched_vessels", [])},
            d.get("target", ""))


def build_rows(patient_dir, pid, wi, study):  # study: "study_a"(PCCT) or "study_b"(EID)
    idx = 0 if study == "study_a" else 1
    rows = []
    for target in ("left", "right"):
        base = os.path.join(patient_dir, "final", target)
        vtxt = os.path.join(base, "common_vessel_sections_volume_summary.txt")
        vjson = os.path.join(base, "common_vessel_sections_summary.json")
        if not (os.path.exists(vtxt) and os.path.exists(vjson)):
            continue
        comp = parse_volume_txt(vtxt)
        lengths, tgt = vessel_lengths(vjson)
        body = "LeftCoronary" if tgt == "left" else "RightCoronary"
        vessels = set()
        for c in comp.values():
            vessels.update(c.keys())
        for v in sorted(vessels):
            lumen = comp.get("lumen", {}).get(v, (0, 0))[idx]
            wall = comp.get("wall", {}).get(v, (0, 0))[idx]      # LumenAndWall region
            calc = comp.get("CALC", {}).get(v, (0, 0))[idx]
            lrnc = comp.get("LRNC", {}).get(v, (0, 0))[idx]
            ncm = comp.get("NonCALCMATX", {}).get(v, (0, 0))[idx]
            rows.append({
                "workitemID": wi, "individualID": pid, "bodySite": body, "location": v,
                "level": "vessel", "Len": lengths.get(v, ""),
                "LumenVol": lumen, "WallVol": wall - lumen, "LumenAndWallVol": wall,
                "CALCVol": calc, "LRNCVol": lrnc, "NonCALCMATXVol": ncm,
                "TotalPlaqueVolume": calc + ncm,
            })
    return rows


def main():
    out_root, pairs_file, out_csv = sys.argv[1], sys.argv[2], sys.argv[3]
    os.makedirs(os.path.join(out_csv, "PCCT"), exist_ok=True)
    os.makedirs(os.path.join(out_csv, "EID"), exist_ok=True)
    n = 0
    for line in open(pairs_file, encoding="utf-8"):
        line = line.strip()
        if not line:
            continue
        pid, pdir, edir = line.split("|")
        patient_dir = os.path.join(out_root, pid)
        if not os.path.isdir(os.path.join(patient_dir, "final")):
            print(f"  SKIP {pid}: no final/ output"); continue
        for study, scan, widir in (("study_a", "PCCT", pdir), ("study_b", "EID", edir)):
            wi = os.path.basename(widir)
            rows = build_rows(patient_dir, pid, wi, study)
            if not rows:
                print(f"  {pid} {scan}: no vessel rows"); continue
            out = os.path.join(out_csv, scan, f"{pid}_workitem_summary_{wi}_subseg.csv")
            with open(out, "w", newline="", encoding="utf-8") as f:
                w = csv.DictWriter(f, fieldnames=HEADER); w.writeheader()
                for r in rows:
                    w.writerow({k: r.get(k, "") for k in HEADER})
        n += 1
        print(f"  {pid}: wrote PCCT+EID ({len(build_rows(patient_dir, pid, os.path.basename(pdir), 'study_a'))} L+R vessels)")
    print(f"Done: {n} patients")


if __name__ == "__main__":
    main()
