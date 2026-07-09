"""Gate 4 v1→v2 shift figures: per-patient Bland-Altman points for v1 (original)
and v2 (2026-07-07), with arrows showing each patient's shift, on BOTH the
untransformed (raw mm³) and log(x+1) scales. One 2-panel figure per endpoint.

Canonical (vessel-overlap) paired data. Writes gate_results/gate4_v1v2_shift/BA_shift_<var>.png.
Run: python scripts/generate_gate4_shift_figs.py
"""
import csv, math, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))
V1 = os.path.join(ROOT, "gate_results_v1_original", "paired_data.csv")
V2 = os.path.join(ROOT, "gate_results", "paired_data.csv")
OUT = os.path.join(ROOT, "gate_results", "gate4_v1v2_shift")

ENDPOINTS = [("LumenVol", "Lumen Volume"), ("WallVol", "Wall Volume"),
             ("VesselVol", "Vessel Volume"), ("CALCVol", "CALC Volume"),
             ("LRNCVol", "LRNC Volume"), ("NonCALCMATXVol", "NonCALC Matrix Volume"),
             ("TotalPlaqueVolume", "Total Plaque Volume")]
C_V1, C_V2, C_ARR = "#93c5fd", "#1d4ed8", "#94a3b8"


def load(path):
    return {r["patient_id"]: r for r in csv.DictReader(open(path, encoding="utf-8"))}


def val(r, k):
    try:
        return float(r[k])
    except Exception:
        return None


def points(rows, var, logt):
    """Return {pid: (mean, diff)} on raw or log(x+1) scale."""
    out = {}
    for pid, r in rows.items():
        p, e = val(r, f"PCCT_{var}"), val(r, f"EID_{var}")
        if p is None or e is None:
            continue
        if logt:
            p, e = math.log(p + 1), math.log(e + 1)
        out[pid] = ((p + e) / 2.0, p - e)
    return out


def panel(ax, d1, d2, scale_label):
    common = sorted(set(d1) & set(d2))
    for pid in common:
        (m1, y1), (m2, y2) = d1[pid], d2[pid]
        ax.annotate("", xy=(m2, y2), xytext=(m1, y1),
                    arrowprops=dict(arrowstyle="->", color=C_ARR, alpha=0.55, lw=0.9))
    x1 = [d1[p][0] for p in common]; y1 = [d1[p][1] for p in common]
    x2 = [d2[p][0] for p in common]; y2 = [d2[p][1] for p in common]
    ax.scatter(x1, y1, c=C_V1, s=45, edgecolors="white", linewidths=0.5, zorder=4, label="v1 (original)")
    ax.scatter(x2, y2, c=C_V2, s=45, edgecolors="white", linewidths=0.5, zorder=5, label="v2 (2026-07-07)")
    b1, b2 = np.mean(y1), np.mean(y2)
    ax.axhline(b1, color=C_V1, ls="--", lw=1.4, label=f"v1 bias {b1:.3f}")
    ax.axhline(b2, color=C_V2, ls="-", lw=1.6, label=f"v2 bias {b2:.3f}")
    ax.axhline(0, color="#64748b", ls=":", lw=0.8)
    ax.set_xlabel(f"Mean of PCCT & EID ({scale_label})", fontsize=9)
    ax.set_ylabel(f"Difference PCCT − EID ({scale_label})", fontsize=9)
    ax.set_title(f"{scale_label}   (N={len(common)}; Δbias {b2 - b1:+.3f})", fontsize=10)
    ax.legend(fontsize=7.5, loc="best", framealpha=0.9)
    ax.grid(True, alpha=0.2)


def main():
    os.makedirs(OUT, exist_ok=True)
    r1, r2 = load(V1), load(V2)
    for var, label in ENDPOINTS:
        fig, (axL, axR) = plt.subplots(1, 2, figsize=(13, 5.4))
        panel(axL, points(r1, var, False), points(r2, var, False), "untransformed, mm³")
        panel(axR, points(r1, var, True), points(r2, var, True), "log(x+1)")
        fig.suptitle(f"Gate 4 — {label}: per-patient v1→v2 shift (canonical)", fontsize=12, y=1.00)
        fig.tight_layout()
        p = os.path.join(OUT, f"BA_shift_{var}.png")
        fig.savefig(p, dpi=150, bbox_inches="tight"); plt.close(fig)
        print("wrote", os.path.relpath(p, ROOT))


if __name__ == "__main__":
    main()
