"""Cohort summary of vessel & plaque characteristics across the paired PCCT/EID scans.

Two tables, one per data vintage (v1 = original, v2 = 2026-07-07). Within each table
PCCT and EID are separate columns (side by side) for the raw (full-vessel, with length
variability) and sub-segment (extent-matched) bases, plus the paired bias Δ = PCCT−EID.
Length is also summarized with its coefficient of variation and the PCCT−EID extent
differential, which the sub-segment removes.

Writes gate_results/cohort_characteristics.csv and .md. Run:
    python scripts/cohort_characteristics.py
"""
import csv
import os
import numpy as np

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))
# vintage -> {basis: paired_data path}
VINTAGES = [
    ("v1 (original)", {
        "raw":     os.path.join(ROOT, "gate_results_v1_original", "paired_data.csv"),
        "sub-seg": os.path.join(ROOT, "gate_results_v1_original", "paired_data_subsegment.csv"),
    }),
    ("v2 (2026-07-07)", {
        "raw":     os.path.join(ROOT, "gate_results", "paired_data.csv"),
        "sub-seg": os.path.join(ROOT, "gate_results", "paired_data_subsegment.csv"),
    }),
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


def summ(p, e):
    """Return (PCCT, EID, Δ) cell strings 'mean±SD (min–max)'; blank if empty."""
    if len(p) == 0:
        return "—", "—", "—"
    d = p - e
    big = p.max() > 20
    nd = 0 if big else 1
    def cell(a, signed=False):
        m = f"{a.mean():+.{nd}f}" if signed else f"{a.mean():.{nd}f}"
        return f"{m}±{a.std(ddof=1):.{nd}f} ({a.min():.{nd}f}–{a.max():.{nd}f})"
    return cell(p), cell(e), cell(d, signed=True)


def pooled_iqr(rows, var):
    """Median [Q1–Q3] over the pooled PCCT+EID values for a metric."""
    p, e = paired(rows, var)
    if len(p) == 0:
        return "—"
    a = np.concatenate([p, e])
    med = np.median(a); q1, q3 = np.percentile(a, [25, 75])
    nd = 0 if a.max() > 20 else 1
    return f"{med:.{nd}f} [{q1:.{nd}f}–{q3:.{nd}f}]"


def main():
    data = {v: {b: load(p) for b, p in bases.items()} for v, bases in VINTAGES}
    ns = {v: {b: len(rows) for b, rows in bases.items()} for v, bases in data.items()}

    # ---- CSV: one row per (vintage, characteristic) ----
    csv_path = os.path.join(ROOT, "gate_results", "cohort_characteristics.csv")
    with open(csv_path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["vintage", "characteristic", "unit",
                    "PCCT_raw", "EID_raw", "delta_raw",
                    "PCCT_subseg", "EID_subseg", "delta_subseg"])
        for vintage, _ in VINTAGES:
            for var, lab, unit in CHARS:
                pr, er, dr = summ(*paired(data[vintage]["raw"], var))
                ps, es, ds = summ(*paired(data[vintage]["sub-seg"], var))
                w.writerow([vintage, lab, unit, pr, er, dr, ps, es, ds])

    # ---- Markdown: one table per vintage, PCCT vs EID side by side ----
    md = ["# Cohort vessel & plaque characteristics — paired PCCT vs EID\n",
          "**Mean ± SD (min–max)** across the paired cohort. **Raw** = full traced vessel (retains "
          "length variability); **sub-seg** = PCCT∩EID extent-matched intersection. Δ = paired PCCT − EID.\n"]
    for vintage, _ in VINTAGES:
        md.append(f"## {vintage}  (raw N={ns[vintage]['raw']}, sub-seg N={ns[vintage]['sub-seg']})\n")
        md.append("| Characteristic | PCCT raw | EID raw | Δ raw | PCCT sub-seg | EID sub-seg | Δ sub-seg |")
        md.append("|---|--:|--:|--:|--:|--:|--:|")
        for var, lab, unit in CHARS:
            pr, er, dr = summ(*paired(data[vintage]["raw"], var))
            ps, es, ds = summ(*paired(data[vintage]["sub-seg"], var))
            md.append(f"| **{lab}** ({unit}) | {pr} | {er} | {dr} | {ps} | {es} | {ds} |")
        md.append("")

    # ---- Dataset characteristics (v2 raw; pooled PCCT+EID, median [IQR]) ----
    v2raw = data["v2 (2026-07-07)"]["raw"]
    n_v2 = ns["v2 (2026-07-07)"]["raw"]
    md.append(f"## Dataset characteristics (v2, raw; N={n_v2}, pooled PCCT+EID, median [Q1–Q3])\n")
    md.append("| Characteristic | Median [Q1–Q3] |")
    md.append("|---|--:|")
    dc_rows = []
    for var, lab, unit in CHARS:
        c = pooled_iqr(v2raw, var)
        md.append(f"| **{lab}** ({unit}) | {c} |")
        dc_rows.append([lab, unit, c])
    md.append("\n*Pooled = PCCT and EID values combined (n = 2 × N per metric). Median [Q1–Q3] "
              "across the pooled distribution.*\n")
    dc_path = os.path.join(ROOT, "gate_results", "dataset_characteristics_pooled.csv")
    with open(dc_path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["characteristic", "unit", "median_IQR_v2raw_pooled"])
        w.writerows(dc_rows)

    # ---- Length variability callout ----
    md.append("## Length variability (the raw-vs-sub-seg point)\n")
    md.append("| Vintage | Basis | PCCT len CV% | EID len CV% | mean \\|PCCT−EID\\| len | % of mean |")
    md.append("|---|---|--:|--:|--:|--:|")
    for vintage, _ in VINTAGES:
        for basis in ("raw", "sub-seg"):
            p, e = paired(data[vintage][basis], "Len")
            if len(p) == 0:
                continue
            d = np.abs(p - e); ml = np.concatenate([p, e]).mean()
            md.append(f"| {vintage} | {basis} | {p.std(ddof=1)/p.mean()*100:.1f}% | "
                      f"{e.std(ddof=1)/e.mean()*100:.1f}% | {d.mean():.0f} mm | {d.mean()/ml*100:.0f}% |")
    md.append("\n*Raw retains a substantial PCCT−EID length differential; the sub-segment "
              "intersection matches extent, so its length differential ≈ 0.*")

    md_path = os.path.join(ROOT, "gate_results", "cohort_characteristics.md")
    open(md_path, "w", encoding="utf-8").write("\n".join(md) + "\n")
    print("\n".join(md))
    print(f"\nWrote {os.path.relpath(csv_path, ROOT)}, {os.path.relpath(dc_path, ROOT)}, "
          f"and {os.path.relpath(md_path, ROOT)}")


if __name__ == "__main__":
    main()
