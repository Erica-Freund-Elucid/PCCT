"""Cohort summary of vessel & plaque characteristics across the paired PCCT/EID scans,
for v1 (original) and v2 (2026-07-07), on the raw (full-vessel, with length variability)
and sub-segment (extent-matched) bases.

For each characteristic reports, across the paired cohort, PCCT mean±SD, EID mean±SD,
and the paired bias Δ = PCCT−EID (mean±SD). Length is included with its coefficient of
variation and the PCCT−EID extent differential, which the sub-segment removes.

Writes gate_results/cohort_characteristics.csv and .md. Run:
    python scripts/cohort_characteristics.py
"""
import csv
import os
import numpy as np

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))
DATASETS = [
    ("v1 raw",      os.path.join(ROOT, "gate_results_v1_original", "paired_data.csv")),
    ("v2 raw",      os.path.join(ROOT, "gate_results", "paired_data.csv")),
    ("v1 sub-seg",  os.path.join(ROOT, "gate_results_v1_original", "paired_data_subsegment.csv")),
    ("v2 sub-seg",  os.path.join(ROOT, "gate_results", "paired_data_subsegment.csv")),
]
# characteristic key -> (label, unit)
CHARS = [
    ("VesselVol",         "Vessel (lumen+wall)", "mm³"),
    ("LumenVol",          "Lumen",               "mm³"),
    ("WallVol",           "Wall",                "mm³"),
    ("TotalPlaqueVolume", "Total plaque",        "mm³"),
    ("CALCVol",           "CALC",                "mm³"),
    ("LRNCVol",           "LRNC",                "mm³"),
    ("NonCALCMATXVol",    "NonCALC matrix",      "mm³"),
    ("Len",               "Vessel length",       "mm"),
]


def load(path):
    if not os.path.isfile(path):
        return []
    return list(csv.DictReader(open(path, encoding="utf-8")))


def paired(rows, var):
    p, e = [], []
    for r in rows:
        try:
            pv, ev = float(r[f"PCCT_{var}"]), float(r[f"EID_{var}"])
        except (KeyError, ValueError, TypeError):
            continue
        p.append(pv); e.append(ev)
    return np.array(p), np.array(e)


def cell(p, e):
    """PCCT / EID / Δ summaries; returns (pcct, eid, delta) mean±SD strings + n."""
    if len(p) == 0:
        return "—", "—", "—", 0
    d = p - e
    fmt = (lambda a: f"{a.mean():.0f}±{a.std(ddof=1):.0f}") if p.max() > 20 else \
          (lambda a: f"{a.mean():.1f}±{a.std(ddof=1):.1f}")
    dd = f"{d.mean():+.0f}±{d.std(ddof=1):.0f}" if abs(d).max() > 20 else f"{d.mean():+.1f}±{d.std(ddof=1):.1f}"
    return fmt(p), fmt(e), dd, len(p)


def main():
    data = {name: load(path) for name, path in DATASETS}
    ns = {name: len(rows) for name, rows in data.items()}

    # CSV: one row per (characteristic, stat)
    csv_path = os.path.join(ROOT, "gate_results", "cohort_characteristics.csv")
    with open(csv_path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["characteristic", "unit", "stat"] + [n for n, _ in DATASETS])
        for var, lab, unit in CHARS:
            cells = {name: cell(*paired(data[name], var)) for name, _ in DATASETS}
            for i, stat in enumerate(("PCCT mean±SD", "EID mean±SD", "Δ PCCT−EID")):
                w.writerow([lab, unit, stat] + [cells[name][i] for name, _ in DATASETS])

    # Markdown table
    md = []
    md.append("# Cohort vessel & plaque characteristics — paired PCCT vs EID\n")
    md.append("Mean ± SD across the paired cohort. **Raw** = full traced vessel "
              "(retains length variability); **sub-seg** = PCCT∩EID extent-matched "
              "intersection. Δ = paired PCCT − EID bias.\n")
    md.append("N: " + ", ".join(f"{n}={ns[n]}" for n, _ in DATASETS) + "\n")
    header = "| Characteristic | Stat | " + " | ".join(n for n, _ in DATASETS) + " |"
    md.append(header)
    md.append("|" + "---|" * (2 + len(DATASETS)))
    for var, lab, unit in CHARS:
        cells = {name: cell(*paired(data[name], var)) for name, _ in DATASETS}
        for i, stat in enumerate(("PCCT", "EID", "Δ (bias)")):
            first = f"**{lab}** ({unit})" if i == 0 else ""
            md.append(f"| {first} | {stat} | " +
                      " | ".join(cells[name][i] for name, _ in DATASETS) + " |")

    # Length variability callout
    md.append("\n## Length variability (the raw-vs-sub-seg point)\n")
    md.append("| Dataset | PCCT len CV% | EID len CV% | mean \\|PCCT−EID\\| len | as % of mean len |")
    md.append("|---|---|---|---|---|")
    for name, _ in DATASETS:
        p, e = paired(data[name], "Len")
        if len(p) == 0:
            continue
        d = np.abs(p - e)
        ml = np.concatenate([p, e]).mean()
        md.append(f"| {name} | {p.std(ddof=1)/p.mean()*100:.1f}% | {e.std(ddof=1)/e.mean()*100:.1f}% | "
                  f"{d.mean():.0f} mm | {d.mean()/ml*100:.0f}% |")
    md.append("\n*Raw retains a substantial PCCT−EID length differential; the sub-segment "
              "intersection matches extent, so its length differential ≈ 0.*")

    md_path = os.path.join(ROOT, "gate_results", "cohort_characteristics.md")
    open(md_path, "w", encoding="utf-8").write("\n".join(md) + "\n")
    print("\n".join(md))
    print(f"\nWrote {os.path.relpath(csv_path, ROOT)} and {os.path.relpath(md_path, ROOT)}")


if __name__ == "__main__":
    main()
