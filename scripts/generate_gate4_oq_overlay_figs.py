"""Gate 4 OQ-overlay Bland-Altman figures — plaque endpoints, on the OQ's basis:
UNTRANSFORMED and LENGTH-NORMALIZED (per-mm). The 730-CVV-040 ATTACHMENT 2
Bland-Altman reference is untransformed + length-normalized, so this is the only
scale on which the OQ bias / bias-CI / LoA can be validly overlaid.

Each figure shows the per-patient PCCT−EID differences (canonical, length-normalized),
the PCCT mean bias + bootstrap 95% CI, and the OQ bias line, OQ bias 95% CI band, and
OQ LoA band — so the CI-overlap acceptance is visible.

Canonical paired data (gate_results/paired_data.csv). Writes
gate_results/bland_altman_plots_oq/BA_oq_<var>.png. Run:
    python scripts/generate_gate4_oq_overlay_figs.py
"""
import csv, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))
# Sub-segment (extent-matched) paired data: on the sub-segment PCCT Len = EID Len, so
# length-normalization is clean and directly comparable to the OQ per-mm reference.
# (Canonical per-mm is confounded by PCCT tracing ~+19% longer than EID.)
PAIRED = os.path.join(ROOT, "gate_results", "paired_data_subsegment.csv")
OUT = os.path.join(ROOT, "gate_results", "bland_altman_plots_oq")
SEED, NBOOT = 12345, 5000

# OQ untransformed, length-normalized BA references (730-CVV-040 ATTACHMENT 2):
# label, bias point, bias 95% CI, LoA (lower, upper)
OQ = {
    "CALCVol":          ("CALC Volume",          0.02,  (-0.01, 0.05), (-0.17, 0.21)),
    "LRNCVol":          ("LRNC Volume",          -0.03, (-0.06, -0.01), (-0.21, 0.14)),
    "NonCALCMATXVol":   ("NonCALC Matrix Volume", -0.07, (-0.19, 0.05), (-0.82, 0.67)),
    "TotalPlaqueVolume": ("Total Plaque Volume",  -0.05, (-0.19, 0.08), (-0.89, 0.78)),
}


def _f(r, k):
    try:
        return float(r[k])
    except Exception:
        return None


def norm_pairs(rows, var):
    """Length-normalized (per-mm) paired values, same filter as the stat report."""
    xs, ys = [], []
    for r in rows:
        p, e = _f(r, f"PCCT_{var}"), _f(r, f"EID_{var}")
        pl, el = _f(r, "PCCT_Len"), _f(r, "EID_Len")
        if None in (p, e, pl, el) or pl <= 0 or el <= 0:
            continue
        xs.append(p / pl); ys.append(e / el)
    return np.array(xs), np.array(ys)


def main():
    os.makedirs(OUT, exist_ok=True)
    rows = list(csv.DictReader(open(PAIRED, encoding="utf-8")))
    rng = np.random.RandomState(SEED)
    for var, (label, oq_bias, oq_ci, oq_loa) in OQ.items():
        pv, ev = norm_pairs(rows, var)
        n = len(pv)
        if n < 2:
            continue
        d = pv - ev
        means = (pv + ev) / 2.0
        bias = float(d.mean())
        bt = [float(np.mean(d[rng.randint(0, n, n)])) for _ in range(NBOOT)]
        b_lo, b_hi = float(np.percentile(bt, 2.5)), float(np.percentile(bt, 97.5))
        overlap = (b_lo <= oq_ci[1]) and (oq_ci[0] <= b_hi)

        fig, ax = plt.subplots(figsize=(7.6, 5.4))
        # OQ bands
        ax.axhspan(oq_loa[0], oq_loa[1], color="#22c55e", alpha=0.07, zorder=0,
                   label="OQ LoA band")
        ax.axhspan(oq_ci[0], oq_ci[1], color="#16a34a", alpha=0.20, zorder=1,
                   label="OQ bias 95% CI")
        ax.axhline(oq_bias, color="#15803d", ls="--", lw=1.6, zorder=3,
                   label=f"OQ bias {oq_bias:+.3f}")
        # PCCT
        ax.scatter(means, d, c="#1d4ed8", s=46, edgecolors="white", linewidths=0.5,
                   zorder=5, label=f"PCCT−EID (N={n})")
        ax.axhline(bias, color="#1d4ed8", ls="-", lw=1.8, zorder=4,
                   label=f"PCCT bias {bias:+.3f} [{b_lo:+.3f},{b_hi:+.3f}]")
        ax.axhspan(b_lo, b_hi, color="#1d4ed8", alpha=0.12, zorder=2)
        ax.axhline(0, color="#64748b", ls=":", lw=0.8)
        ax.set_xlabel("Mean of PCCT & EID  (length-normalized, per-mm)", fontsize=9)
        ax.set_ylabel("Difference PCCT − EID  (per-mm)", fontsize=9)
        verdict = "OVERLAP (pass)" if overlap else "NO OVERLAP (bias exceeds OQ)"
        ax.set_title(f"{label} — sub-segment (extent-matched), length-normalized vs OQ\n{verdict}",
                     fontsize=11)
        ax.legend(fontsize=7.5, loc="best", framealpha=0.92)
        ax.grid(True, alpha=0.2)
        fig.tight_layout()
        p = os.path.join(OUT, f"BA_oq_{var}.png")
        fig.savefig(p, dpi=150, bbox_inches="tight"); plt.close(fig)
        print("wrote", os.path.relpath(p, ROOT), "overlap=" + str(overlap))


if __name__ == "__main__":
    main()
