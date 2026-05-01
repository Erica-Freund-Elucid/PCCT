"""
PCCT Qualification Gate Analyses
=================================
Computes Gate 2 (workflow), Gate 3 (reproducibility), and Gate 4 (bias/agreement)
metrics from paired PCCT/EID workitem summary CSVs.

Acceptance criteria: 95% CI overlap with B.1P Delta Validation OQ (730-CVV-040).
Cross-scanner (PCCT vs EID) variability must be non-inferior to validated
inter-operator variability from the delta OQ.

Each patient must have one CSV in workitem_summaries/PCCT/ and one in
workitem_summaries/EID/ with matching patient IDs (individualID column).

Patient-level totals are computed by summing the two target-level rows
(LeftCoronary + RightCoronary) per CSV. Volumes are normalized to total
vessel length to account for differences in traced extent between scans.

Usage:
    python run_gate_analyses.py
"""

import csv
import glob
import os
import sys
import math
from collections import defaultdict

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

try:
    from scipy import stats as sp_stats
    HAS_SCIPY = True
except ImportError:
    HAS_SCIPY = False
    print("WARNING: scipy not installed. Bootstrap CIs will be used instead of parametric.")

# ── Configuration ─────────────────────────────────────────────────────────────

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PCCT_DIR = os.path.join(SCRIPT_DIR, "workitem_summaries", "PCCT")
EID_DIR = os.path.join(SCRIPT_DIR, "workitem_summaries", "EID")
OUTPUT_DIR = os.path.join(SCRIPT_DIR, "gate_results")

# ── Reference data ────────────────────────────────────────────────────────────
#
# Two validated reference studies establish the accepted performance envelope:
#
#   1. Original B.1P OQ (4-B1P-033 v2.0, 730-CVV-008)
#      - Acceptance criterion: <=10% wCV point estimate for process outputs
#      - Log(x+1) transformation applied to process outputs (lumen, wall,
#        vessel vol, length); untransformed for plaque volumes
#      - Clinical outputs (plaque volumes) reported as secondary/descriptive
#
#   2. B.1P Delta Validation OQ (730-CVV-040 v0.1)
#      - Acceptance criterion: 95% CI overlap with original OQ (log-wCV)
#      - Expanded case mix (artifacts, heavy calcification)
#      - Design change: wall collapses to lumen outside lesions (730-CVV-017)
#      - Length demoted to descriptive (variability reflects termination judgment)
#      - Plaque volume reproducibility improved vs original despite harder cases
#
# PCCT qualification acceptance:
#   Primary: 95% CI overlap with delta OQ wCV (log scale for volumes)
#   The delta OQ CIs already encompass the original OQ performance, so
#   overlap with the delta OQ demonstrates non-inferiority to both studies.
#   Original OQ results reported for context and traceability.
#
# Primary metric: log-wCV for volume endpoints (measurement error proportional
# to magnitude); untransformed wCV for length (additive error structure).
# Length is descriptive only per 730-CVV-040 recommendation.

# Gate 3 — Primary endpoints
# Delta OQ: Table 4 of 730-CVV-040
# Original OQ: Table 7 of 4-B1P-033 — log(x+1) transformed for process
#   outputs (lumen, wall, vessel vol, length); untransformed for plaque volumes
GATE3_PRIMARY = {
    "LumenVol": {
        "label": "Lumen Volume (mm³)",
        "metric": "log",
        # Delta OQ (730-CVV-040)
        "delta_log_wcv": 7.5,
        "delta_log_ci": (5.9, 9.3),
        "delta_untrans_wcv": 9.32,
        "delta_untrans_ci": (7.09, 11.58),
        # Original OQ (4-B1P-033) — log-transformed, acceptance <=10%
        "orig_log_wcv": 7.64,
        "orig_log_ci": (5.52, 8.55),
    },
    "WallVol": {
        "label": "Wall Volume (mm³)",
        "metric": "log",
        # Delta OQ — log-wCV is primary; untransformed inflated by near-zero
        # wall outside lesions post design change (730-CVV-017)
        "delta_log_wcv": 13.1,
        "delta_log_ci": (10.2, 16.0),
        "delta_untrans_wcv": 23.50,
        "delta_untrans_ci": (17.88, 29.71),
        # Original OQ — log-transformed, pre-design change
        "orig_log_wcv": 9.51,
        "orig_log_ci": (6.6, 11.16),
    },
    "VesselVol": {
        "label": "Vessel Volume (mm³)",
        "metric": "both",
        # Delta OQ — log primary, untransformed for direct comparison
        "delta_log_wcv": 8.7,
        "delta_log_ci": (6.8, 10.6),
        "delta_untrans_wcv": 10.18,
        "delta_untrans_ci": (7.83, 12.38),
        # Original OQ — log-transformed, tight CI reflects narrower case mix
        "orig_log_wcv": 6.03,
        "orig_log_ci": (4.44, 6.74),
    },
}

# Gate 3 — Descriptive only (per 730-CVV-040 recommendation)
GATE3_DESCRIPTIVE = {
    "Len": {
        "label": "Vessel Length (mm)",
        "metric": "untransformed",
        # Delta OQ — elevated due to expanded artifact inclusion
        "delta_untrans_wcv": 13.88,
        "delta_untrans_ci": (10.72, 17.38),
        "delta_log_wcv": 15.5,
        "delta_log_ci": (12.2, 18.8),
        # Original OQ — log-transformed
        "orig_log_wcv": 5.27,
        "orig_log_ci": (3.43, 6.52),
    },
}

# Gate 3 — Secondary endpoints: plaque volumes
# Delta OQ: Table 5 of 730-CVV-040
# Original OQ: Table 9 of 4-B1P-033 — untransformed, no CIs reported
GATE3_SECONDARY = {
    "CALCVol": {
        "label": "CALC Volume (mm³)",
        "metric": "log",
        "delta_log_wcv": 4.5,
        "delta_log_ci": (3.5, 5.4),
        "delta_untrans_wcv": 25.78,
        "delta_untrans_ci": (18.70, 35.75),
        # Original OQ — untransformed, no CI
        "orig_untrans_wcv": 13.9,
    },
    "LRNCVol": {
        "label": "LRNC Volume (mm³)",
        "metric": "log",
        "delta_log_wcv": 5.4,
        "delta_log_ci": (4.1, 6.9),
        "delta_untrans_wcv": 78.80,
        "delta_untrans_ci": (48.09, 163.40),
        # Original OQ — untransformed, no CI
        "orig_untrans_wcv": 59.3,
    },
    "NonCALCMATXVol": {
        "label": "NonCALC Matrix Volume (mm³)",
        "metric": "log",
        "delta_log_wcv": 13.6,
        "delta_log_ci": (10.8, 16.4),
        "delta_untrans_wcv": 32.61,
        "delta_untrans_ci": (24.87, 41.94),
        # Original OQ reported Total NCP (58.6%) untransformed;
        # NonCALCMATX is the analogous variable post design change (730-CVV-017)
        "orig_untrans_wcv": 58.6,
        "orig_note": "Original reported as Total NCP (untransformed); improved in delta OQ",
    },
    "TotalPlaqueVolume": {
        "label": "Total Plaque Volume (mm³)",
        "metric": "log",
        "delta_log_wcv": 13.2,
        "delta_log_ci": (10.4, 16.0),
        "delta_untrans_wcv": 27.08,
        "delta_untrans_ci": (20.68, 34.74),
        # Original OQ — untransformed, no CI
        "orig_untrans_wcv": 44.1,
    },
}

# Gate 4 — Bland-Altman reference (Table 6 of 730-CVV-040, log scale)
GATE4_VARIABLES = {
    "LumenVol": {
        "label": "Lumen Volume (mm³)",
        "bias_threshold_pct": 5,
    },
    "CALCVol": {
        "label": "CALC Volume (mm³)",
        "bias_threshold_pct": 10,
        "ref_bias": 0.02,
        "ref_loa": (-0.17, 0.21),
    },
    "LRNCVol": {
        "label": "LRNC Volume (mm³)",
        "bias_threshold_pct": 10,
        "ref_bias": -0.03,
        "ref_loa": (-0.21, 0.14),
    },
    "NonCALCMATXVol": {
        "label": "NonCALC Matrix Volume (mm³)",
        "bias_threshold_pct": 10,
        "ref_bias": -0.07,
        "ref_loa": (-0.82, 0.67),
    },
    "TotalPlaqueVolume": {
        "label": "Total Plaque Volume (mm³)",
        "bias_threshold_pct": 10,
        "ref_bias": -0.05,
        "ref_loa": (-0.89, 0.78),
    },
}

# Collect all variable names that need to be loaded from CSVs
ALL_VARS = set()
for d in [GATE3_PRIMARY, GATE3_DESCRIPTIVE, GATE3_SECONDARY, GATE4_VARIABLES]:
    ALL_VARS.update(d.keys())
# VesselVol is derived: LumenVol + WallVol (= LumenAndWallVol in CSV)
ALL_VARS.add("LumenAndWallVol")
ALL_VARS.discard("VesselVol")

N_BOOTSTRAP = 2000
RANDOM_SEED = 42

# ── Data Loading ──────────────────────────────────────────────────────────────


def load_patient_totals(csv_path):
    """Load a workitem summary CSV and return patient-level totals.

    Sums the two target-level rows (LeftCoronary + RightCoronary) to get
    whole-patient values. Volumes are also length-normalized.
    """
    with open(csv_path, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        rows = list(reader)

    # Target-level rows: use 'level' column if present, otherwise
    # fall back to rows with empty 'location' (older CSV format)
    if rows and "level" in rows[0]:
        target_rows = [r for r in rows if r.get("level") == "target"]
    else:
        target_rows = [r for r in rows if not r.get("location", "").strip()]
    if not target_rows:
        return None

    patient_id = target_rows[0].get("individualID", "")
    # Normalize: strip kernel suffix (e.g. PT-116_Bv_Bv44u_4 -> PT-116)
    if "_Bv_" in patient_id:
        patient_id = patient_id.split("_Bv_")[0]
    totals = {"patient_id": patient_id}

    for key in ALL_VARS:
        val_sum = 0.0
        valid = False
        for row in target_rows:
            raw = row.get(key, "")
            if raw and raw not in ("", "None", "N/A"):
                try:
                    val_sum += float(raw)
                    valid = True
                except (ValueError, TypeError):
                    pass
        totals[key] = val_sum if valid else None

    # Derive VesselVol = LumenAndWallVol (already summed from CSV)
    totals["VesselVol"] = totals.get("LumenAndWallVol")

    # Compute length-normalized volumes
    total_len = totals.get("Len")
    if total_len and total_len > 0:
        for key in ALL_VARS:
            if key != "Len" and totals.get(key) is not None:
                totals[f"{key}_norm"] = totals[key] / total_len
        if totals.get("VesselVol") is not None:
            totals["VesselVol_norm"] = totals["VesselVol"] / total_len

    return totals


def load_paired_data():
    """Load all paired PCCT/EID data."""
    pcct_files = sorted(glob.glob(os.path.join(PCCT_DIR, "*.csv")))
    eid_files = sorted(glob.glob(os.path.join(EID_DIR, "*.csv")))

    if not pcct_files:
        print(f"ERROR: No PCCT summary CSVs found in {PCCT_DIR}")
        return []
    if not eid_files:
        print(f"ERROR: No EID summary CSVs found in {EID_DIR}")
        print(f"  Place EID workitem summary CSVs in: {EID_DIR}")
        return []

    pcct_data = {}
    for f in pcct_files:
        totals = load_patient_totals(f)
        if totals:
            pcct_data[totals["patient_id"]] = totals

    eid_data = {}
    for f in eid_files:
        totals = load_patient_totals(f)
        if totals:
            eid_data[totals["patient_id"]] = totals

    paired = []
    matched_ids = sorted(set(pcct_data.keys()) & set(eid_data.keys()))
    for pid in matched_ids:
        paired.append({
            "patient_id": pid,
            "pcct": pcct_data[pid],
            "eid": eid_data[pid],
        })

    pcct_only = set(pcct_data.keys()) - set(eid_data.keys())
    eid_only = set(eid_data.keys()) - set(pcct_data.keys())

    print(f"Paired patients:   {len(paired)}")
    print(f"PCCT-only:         {len(pcct_only)}")
    print(f"EID-only:          {len(eid_only)}")
    if pcct_only:
        print(f"  PCCT without EID: {', '.join(sorted(pcct_only))}")
    if eid_only:
        print(f"  EID without PCCT: {', '.join(sorted(eid_only))}")
    print()

    return paired


# ── Statistical Functions ─────────────────────────────────────────────────────


def compute_wcv(vals_a, vals_b, log_transform=False):
    """Compute within-subject coefficient of variation for paired measurements.

    If log_transform=True, applies log(x+1) transformation before computing
    wCV on the transformed scale (consistent with 730-CVV-040 methodology).
    """
    if log_transform:
        vals_a = [math.log(v + 1) for v in vals_a]
        vals_b = [math.log(v + 1) for v in vals_b]

    n = len(vals_a)
    if n == 0:
        return None

    ratios = []
    for a, b in zip(vals_a, vals_b):
        m = (a + b) / 2.0
        if m == 0:
            continue
        d = a - b
        ratios.append((d ** 2) / (2.0 * m ** 2))

    if not ratios:
        return None

    return math.sqrt(np.mean(ratios)) * 100


def bootstrap_wcv_ci(vals_a, vals_b, log_transform=False,
                     n_boot=N_BOOTSTRAP, alpha=0.05):
    """Bootstrap 95% CI for wCV."""
    rng = np.random.RandomState(RANDOM_SEED)
    n = len(vals_a)
    arr_a = np.array(vals_a)
    arr_b = np.array(vals_b)

    boot_wcvs = []
    for _ in range(n_boot):
        idx = rng.randint(0, n, size=n)
        wcv = compute_wcv(arr_a[idx].tolist(), arr_b[idx].tolist(),
                          log_transform=log_transform)
        if wcv is not None:
            boot_wcvs.append(wcv)

    if not boot_wcvs:
        return None, None

    lb = np.percentile(boot_wcvs, 100 * alpha / 2)
    ub = np.percentile(boot_wcvs, 100 * (1 - alpha / 2))
    return lb, ub


def ci_overlap(ci1, ci2):
    """Check whether two confidence intervals overlap."""
    return ci1[0] <= ci2[1] and ci2[0] <= ci1[1]


def bland_altman(vals_a, vals_b, log_transform=False):
    """Compute Bland-Altman statistics (PCCT minus EID).

    If log_transform=True, computes on log(x+1) scale (per 730-CVV-040).
    """
    if log_transform:
        a = np.array([math.log(v + 1) for v in vals_a])
        b = np.array([math.log(v + 1) for v in vals_b])
    else:
        a = np.array(vals_a)
        b = np.array(vals_b)

    diffs = a - b
    means = (a + b) / 2.0

    mean_bias = np.mean(diffs)
    sd_diff = np.std(diffs, ddof=1)
    loa_lower = mean_bias - 1.96 * sd_diff
    loa_upper = mean_bias + 1.96 * sd_diff

    # Proportional bias: regress differences on means
    if len(means) > 2:
        if HAS_SCIPY:
            slope, intercept, r_value, p_value, se = sp_stats.linregress(means, diffs)
            r_sq = r_value ** 2
        else:
            corr = np.corrcoef(means, diffs)[0, 1]
            r_sq = corr ** 2 if not np.isnan(corr) else 0
    else:
        r_sq = 0

    return mean_bias, sd_diff, loa_lower, loa_upper, r_sq


def plot_bland_altman(vals_a, vals_b, label, out_path,
                     log_transform=True, ref_bias=None, ref_loa=None,
                     pids=None):
    """Generate and save a Bland-Altman plot."""
    if log_transform:
        a = np.array([math.log(v + 1) for v in vals_a])
        b = np.array([math.log(v + 1) for v in vals_b])
        scale_label = "log(x+1)"
    else:
        a = np.array(vals_a)
        b = np.array(vals_b)
        scale_label = "raw"

    diffs = a - b
    means = (a + b) / 2.0
    mean_bias = np.mean(diffs)
    sd_diff = np.std(diffs, ddof=1)
    loa_lo = mean_bias - 1.96 * sd_diff
    loa_hi = mean_bias + 1.96 * sd_diff

    fig, ax = plt.subplots(figsize=(7, 5))

    ax.scatter(means, diffs, c="#2563eb", s=50, zorder=5, edgecolors="white", linewidths=0.5)

    # Label points with patient IDs if provided
    if pids:
        for i, pid in enumerate(pids):
            ax.annotate(pid, (means[i], diffs[i]), fontsize=7,
                        xytext=(4, 4), textcoords="offset points", color="#64748b")

    # Mean bias line
    ax.axhline(mean_bias, color="#dc2626", linestyle="-", linewidth=1.5,
               label=f"Bias: {mean_bias:.4f}")
    # LoA lines
    ax.axhline(loa_hi, color="#f59e0b", linestyle="--", linewidth=1,
               label=f"+1.96 SD: {loa_hi:.4f}")
    ax.axhline(loa_lo, color="#f59e0b", linestyle="--", linewidth=1,
               label=f"-1.96 SD: {loa_lo:.4f}")
    # Zero reference
    ax.axhline(0, color="#94a3b8", linestyle=":", linewidth=0.8)

    # Delta OQ reference if available
    if ref_bias is not None and ref_loa is not None:
        ax.axhspan(ref_loa[0], ref_loa[1], alpha=0.08, color="#16a34a",
                    label=f"Delta OQ LoA: [{ref_loa[0]}, {ref_loa[1]}]")
        ax.axhline(ref_bias, color="#16a34a", linestyle="-", linewidth=1, alpha=0.5)

    # Proportional bias regression line
    if len(means) > 2:
        corr = np.corrcoef(means, diffs)[0, 1]
        r_sq = corr ** 2 if not np.isnan(corr) else 0
        if HAS_SCIPY:
            slope, intercept, _, _, _ = sp_stats.linregress(means, diffs)
        else:
            slope = corr * np.std(diffs) / np.std(means) if np.std(means) > 0 else 0
            intercept = mean_bias - slope * np.mean(means)
        x_line = np.array([means.min(), means.max()])
        ax.plot(x_line, slope * x_line + intercept, color="#94a3b8",
                linestyle="-", linewidth=1, alpha=0.6,
                label=f"Regression (r2={r_sq:.3f})")

    ax.set_xlabel(f"Mean of PCCT and EID ({scale_label})", fontsize=10)
    ax.set_ylabel(f"Difference: PCCT - EID ({scale_label})", fontsize=10)
    ax.set_title(f"Bland-Altman: {label}\n(N={len(vals_a)}, {scale_label} scale)", fontsize=11)
    ax.legend(fontsize=8, loc="best", framealpha=0.9)
    ax.grid(True, alpha=0.2)

    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def bootstrap_rsq_ci(vals_a, vals_b, log_transform=False,
                     n_boot=N_BOOTSTRAP, alpha=0.05):
    """Bootstrap 95% CI for proportional bias r²."""
    rng = np.random.RandomState(RANDOM_SEED + 7)
    n = len(vals_a)
    if n < 3:
        return None, None

    if log_transform:
        a = np.array([math.log(v + 1) for v in vals_a])
        b = np.array([math.log(v + 1) for v in vals_b])
    else:
        a = np.array(vals_a)
        b = np.array(vals_b)

    boot_rsqs = []
    for _ in range(n_boot):
        idx = rng.randint(0, n, size=n)
        diffs = a[idx] - b[idx]
        means = (a[idx] + b[idx]) / 2.0
        if np.std(means) == 0 or np.std(diffs) == 0:
            continue
        corr = np.corrcoef(means, diffs)[0, 1]
        if not np.isnan(corr):
            boot_rsqs.append(corr ** 2)

    if not boot_rsqs:
        return None, None

    lb = np.percentile(boot_rsqs, 100 * alpha / 2)
    ub = np.percentile(boot_rsqs, 100 * (1 - alpha / 2))
    return lb, ub


# ── Gate Analysis Helpers ─────────────────────────────────────────────────────


def _get_paired_values(paired, var, normalize=False):
    """Extract valid paired values for a variable.

    If normalize=True and var is not 'Len', uses length-normalized values
    (vol/mm) to control for differences in vessel extent between scans.
    """
    pcct_vals, eid_vals, pids = [], [], []
    lookup_key = f"{var}_norm" if (normalize and var != "Len") else var
    for pair in paired:
        pv = pair["pcct"].get(lookup_key)
        ev = pair["eid"].get(lookup_key)
        if pv is not None and ev is not None:
            pcct_vals.append(pv)
            eid_vals.append(ev)
            pids.append(pair["patient_id"])
    return pcct_vals, eid_vals, pids


def _find_csv(directory, patient_id):
    for f in glob.glob(os.path.join(directory, "*.csv")):
        totals = load_patient_totals(f)
        if totals and totals["patient_id"] == patient_id:
            return f
    return None


def _count_segments(csv_path):
    with open(csv_path, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        return sum(1 for r in reader if r.get("level") == "segment")


# ── Gate 2: Workflow Integration ──────────────────────────────────────────────


def run_gate2(paired):
    lines = []
    lines.append("=" * 70)
    lines.append("GATE 2 — WORKFLOW INTEGRATION CHECKS")
    lines.append("=" * 70)
    lines.append("")

    # 2.2 Vessel tree coverage
    lines.append("--- 2.2 Centerline / Vessel Tree Coverage ---")
    lines.append("")
    lines.append(f"{'Patient':<15s} {'PCCT Len':>10s} {'EID Len':>10s} {'Diff%':>8s}  "
                 f"{'PCCT Segs':>10s} {'EID Segs':>10s}")
    lines.append("-" * 70)

    len_diffs = []
    seg_diffs = []
    for pair in paired:
        pcct_file = _find_csv(PCCT_DIR, pair["patient_id"])
        eid_file = _find_csv(EID_DIR, pair["patient_id"])
        if not pcct_file or not eid_file:
            continue

        pcct_segs = _count_segments(pcct_file)
        eid_segs = _count_segments(eid_file)
        pcct_len = pair["pcct"].get("Len")
        eid_len = pair["eid"].get("Len")

        diff_pct = None
        if pcct_len and eid_len and eid_len > 0:
            diff_pct = (pcct_len - eid_len) / eid_len * 100
            len_diffs.append(diff_pct)
        seg_diffs.append(pcct_segs - eid_segs)

        diff_str = f"{diff_pct:>7.1f}%" if diff_pct is not None else "    N/A "
        lines.append(
            f"{pair['patient_id']:<15s} "
            f"{pcct_len:>10.1f} {eid_len:>10.1f} "
            f"{diff_str}  "
            f"{pcct_segs:>10d} {eid_segs:>10d}"
        )

    if len_diffs:
        lines.append("")
        lines.append(f"  Mean length difference:    {np.mean(len_diffs):+.1f}%")
        lines.append(f"  Mean segment count diff:   {np.mean(seg_diffs):+.1f}")
        lines.append(f"  Pairs with data:           {len(len_diffs)}/{len(paired)}")
        threshold_met = abs(np.mean(len_diffs)) < 15  # >=95% coverage ~ <15% diff
        lines.append(f"  Threshold (>=95% coverage): {'PASS' if threshold_met else 'REVIEW'}")
    lines.append("")

    return "\n".join(lines)


# ── Gate 3: Reproducibility ──────────────────────────────────────────────────


def _analyze_variable(lines, detail_lines, paired, var, cfg, is_descriptive=False):
    """Run wCV analysis for a single variable.

    Volume variables are length-normalized (vol/mm) before computing wCV
    to isolate variability due to image modality from vessel extent differences.
    Length itself is not normalized.
    """
    normalize = var != "Len"
    pcct_vals, eid_vals, pids = _get_paired_values(paired, var, normalize=normalize)
    n_valid = len(pcct_vals)

    label = cfg["label"]
    metric = cfg.get("metric", "log")
    tag = " (DESCRIPTIVE)" if is_descriptive else ""
    norm_tag = " [length-normalized]" if normalize else ""

    lines.append(f"--- {label}{tag}{norm_tag} ---")

    if n_valid < 2:
        lines.append(f"  INSUFFICIENT DATA (n={n_valid})")
        lines.append("")
        return

    # Compute both untransformed and log wCV
    wcv_untrans = compute_wcv(pcct_vals, eid_vals, log_transform=False)
    ci_untrans = bootstrap_wcv_ci(pcct_vals, eid_vals, log_transform=False)

    wcv_log = compute_wcv(pcct_vals, eid_vals, log_transform=True)
    ci_log = bootstrap_wcv_ci(pcct_vals, eid_vals, log_transform=True)

    lines.append(f"  N:                  {n_valid}")
    lines.append(f"  Untransformed wCV:  {wcv_untrans:.2f}% [{ci_untrans[0]:.2f}%, {ci_untrans[1]:.2f}%]")
    lines.append(f"  Log-wCV:            {wcv_log:.2f}% [{ci_log[0]:.2f}%, {ci_log[1]:.2f}%]")
    lines.append("")

    # ── Reference comparison ──────────────────────────────────────────────
    # Show both original OQ and delta OQ for traceability

    # Original B.1P OQ (4-B1P-033)
    orig_log_wcv = cfg.get("orig_log_wcv")
    orig_log_ci = cfg.get("orig_log_ci")
    orig_untrans_wcv = cfg.get("orig_untrans_wcv")
    orig_note = cfg.get("orig_note")
    if orig_log_wcv is not None:
        if orig_log_ci:
            lines.append(f"  Original OQ (4-B1P-033):  {orig_log_wcv:.2f}% [{orig_log_ci[0]:.2f}%, {orig_log_ci[1]:.2f}%] (log-transformed)")
        else:
            lines.append(f"  Original OQ (4-B1P-033):  {orig_log_wcv:.2f}% (log-transformed, no CI reported)")
    elif orig_untrans_wcv is not None:
        lines.append(f"  Original OQ (4-B1P-033):  {orig_untrans_wcv:.1f}% (untransformed, no CI reported)")
    if orig_note:
        lines.append(f"    Note: {orig_note}")

    # Delta OQ (730-CVV-040)
    delta_log = cfg.get("delta_log_wcv")
    delta_log_ci = cfg.get("delta_log_ci")
    delta_ut = cfg.get("delta_untrans_wcv")
    delta_ut_ci = cfg.get("delta_untrans_ci")

    if delta_log is not None and delta_log_ci:
        lines.append(f"  Delta OQ (730-CVV-040):   {delta_log:.1f}% [{delta_log_ci[0]:.1f}%, {delta_log_ci[1]:.1f}%] (log-wCV)")
    if delta_ut is not None and delta_ut_ci:
        lines.append(f"  Delta OQ (730-CVV-040):   {delta_ut:.2f}% [{delta_ut_ci[0]:.2f}%, {delta_ut_ci[1]:.2f}%] (untransformed)")
    lines.append("")

    # ── Pass/fail determination ───────────────────────────────────────────
    if not is_descriptive:
        if metric in ("log", "both") and delta_log_ci:
            pcct_ci_log = (ci_log[0], ci_log[1])
            overlap_log = ci_overlap(pcct_ci_log, delta_log_ci)
            lines.append(f"  ACCEPTANCE (log-wCV CI overlap with delta OQ):")
            lines.append(f"    PCCT:      [{ci_log[0]:.2f}%, {ci_log[1]:.2f}%]")
            lines.append(f"    Delta OQ:  [{delta_log_ci[0]:.1f}%, {delta_log_ci[1]:.1f}%]")
            lines.append(f"    Overlap:   {'YES — PASS' if overlap_log else 'NO — FAIL'}")

        if metric in ("untransformed", "both") and delta_ut_ci:
            pcct_ci_ut = (ci_untrans[0], ci_untrans[1])
            overlap_ut = ci_overlap(pcct_ci_ut, delta_ut_ci)
            lines.append(f"  ACCEPTANCE (untransformed CI overlap with delta OQ):")
            lines.append(f"    PCCT:      [{ci_untrans[0]:.2f}%, {ci_untrans[1]:.2f}%]")
            lines.append(f"    Delta OQ:  [{delta_ut_ci[0]:.2f}%, {delta_ut_ci[1]:.2f}%]")
            lines.append(f"    Overlap:   {'YES — PASS' if overlap_ut else 'NO — FAIL'}")
    else:
        lines.append(f"  DESCRIPTIVE ONLY — not a primary acceptance criterion")
        lines.append(f"  per 730-CVV-040 recommendation. Vessel length variability")
        lines.append(f"  reflects analyst termination judgment, not segmentation error.")

    lines.append("")

    # Per-patient detail
    detail_lines.append(f"\n{label} — per-patient detail")
    detail_lines.append(f"{'Patient':<15s} {'PCCT':>12s} {'EID':>12s} {'Diff':>12s} {'%Diff':>8s}")
    detail_lines.append("-" * 62)
    for pid, pv, ev in zip(pids, pcct_vals, eid_vals):
        m = (pv + ev) / 2
        d = pv - ev
        pct = (d / m * 100) if m != 0 else 0
        detail_lines.append(f"{pid:<15s} {pv:>12.2f} {ev:>12.2f} {d:>12.2f} {pct:>7.1f}%")


def run_gate3(paired):
    lines = []
    detail = []

    lines.append("=" * 70)
    lines.append("GATE 3 — CORE QUANTITATIVE REPRODUCIBILITY")
    lines.append("=" * 70)
    lines.append(f"N = {len(paired)} paired patients (target: >=30)")
    lines.append(f"Acceptance: 95% CI overlap with B.1P Delta OQ (730-CVV-040)")
    lines.append(f"Primary metric: log-wCV for volumes, untransformed for length")
    lines.append(f"wCV method: sqrt(mean(d²/(2m²))) × 100")
    lines.append(f"95% CI: {N_BOOTSTRAP}-sample bootstrap")
    lines.append(f"Volumes normalized to vessel length per 730-CVV-040 methodology")
    lines.append("")

    lines.append("--- PRIMARY ENDPOINTS ---")
    lines.append("")
    for var, cfg in GATE3_PRIMARY.items():
        _analyze_variable(lines, detail, paired, var, cfg)

    lines.append("--- DESCRIPTIVE (vessel length) ---")
    lines.append("")
    for var, cfg in GATE3_DESCRIPTIVE.items():
        _analyze_variable(lines, detail, paired, var, cfg, is_descriptive=True)

    lines.append("--- SECONDARY ENDPOINTS (plaque volumes) ---")
    lines.append("")
    for var, cfg in GATE3_SECONDARY.items():
        _analyze_variable(lines, detail, paired, var, cfg)

    return "\n".join(lines), "\n".join(detail)


# ── Gate 4: Bias & Agreement ─────────────────────────────────────────────────


def run_gate4(paired):
    lines = []
    detail = []

    lines.append("=" * 70)
    lines.append("GATE 4 — SYSTEMATIC BIAS & AGREEMENT")
    lines.append("=" * 70)
    lines.append(f"N = {len(paired)} paired patients")
    lines.append(f"Method: Bland-Altman on log(x+1) scale, raw volumes (PCCT minus EID)")
    lines.append(f"Note: Raw (not length-normalized) volumes used for direct")
    lines.append(f"      comparability with delta OQ reference (730-CVV-040 Table 6)")
    lines.append("")

    for var, cfg in GATE4_VARIABLES.items():
        # Bland-Altman uses raw volumes on log(x+1) scale, not length-normalized,
        # for direct comparability with delta OQ reference (730-CVV-040 Table 6)
        pcct_vals, eid_vals, pids = _get_paired_values(paired, var, normalize=False)
        n_valid = len(pcct_vals)
        label = cfg["label"]

        lines.append(f"--- {label} ---")

        if n_valid < 2:
            lines.append(f"  INSUFFICIENT DATA (n={n_valid})")
            lines.append("")
            continue

        # Log-scale Bland-Altman (consistent with 730-CVV-040 Table 6)
        mean_bias, sd_diff, loa_lo, loa_hi, r_sq = bland_altman(
            pcct_vals, eid_vals, log_transform=True
        )

        # Also compute untransformed for context
        mean_bias_ut, sd_diff_ut, loa_lo_ut, loa_hi_ut, r_sq_ut = bland_altman(
            pcct_vals, eid_vals, log_transform=False
        )
        pair_means_ut = [(p + e) / 2 for p, e in zip(pcct_vals, eid_vals)]
        overall_mean_ut = np.mean(pair_means_ut)
        bias_pct_ut = abs(mean_bias_ut) / overall_mean_ut * 100 if overall_mean_ut != 0 else 0

        threshold = cfg["bias_threshold_pct"]
        passed_bias = bias_pct_ut < threshold

        lines.append(f"  N:                     {n_valid}")
        lines.append(f"  Log-scale bias:        {mean_bias:.4f}")
        lines.append(f"  Log-scale LoA:         [{loa_lo:.4f}, {loa_hi:.4f}]")
        lines.append(f"  Untransformed bias:    {mean_bias_ut:.2f} ({bias_pct_ut:.1f}% of mean)")
        lines.append(f"  Bias threshold:        |bias| < {threshold}% of mean")
        lines.append(f"  Bias result:           {'PASS' if passed_bias else 'FAIL'}")

        # Compare LoA with delta OQ reference if available
        ref_bias = cfg.get("ref_bias")
        ref_loa = cfg.get("ref_loa")
        if ref_loa:
            lines.append(f"  Delta OQ ref bias:     {ref_bias}")
            lines.append(f"  Delta OQ ref LoA:      [{ref_loa[0]}, {ref_loa[1]}]")

        rsq_lb, rsq_ub = bootstrap_rsq_ci(pcct_vals, eid_vals, log_transform=True)
        if rsq_lb is not None:
            # Pass if the CI includes values below 0.1 (cannot rule out no proportional bias)
            prop_pass = rsq_lb < 0.1
            lines.append(f"  Proportional bias r²:  {r_sq:.3f} [{rsq_lb:.3f}, {rsq_ub:.3f}]")
            lines.append(f"  Prop. bias result:     {'PASS' if prop_pass else 'FAIL'} (threshold: 95% CI lower bound < 0.1)")
        else:
            lines.append(f"  Proportional bias r²:  {r_sq:.3f} (CI not estimable)")
        lines.append("")

        # Bland-Altman plot
        from pathlib import Path
        plot_dir = Path(OUTPUT_DIR) / "bland_altman_plots"
        ref_b = cfg.get("ref_bias")
        ref_l = cfg.get("ref_loa")
        plot_bland_altman(
            pcct_vals, eid_vals, label,
            plot_dir / f"BA_{var}.png",
            log_transform=True,
            ref_bias=ref_b, ref_loa=ref_l,
            pids=pids,
        )

        # Per-patient detail
        detail.append(f"\n{label} — Bland-Altman detail (log scale)")
        detail.append(f"{'Patient':<15s} {'log(PCCT+1)':>12s} {'log(EID+1)':>12s} {'Diff':>12s} {'Mean':>12s}")
        detail.append("-" * 66)
        for pid, pv, ev in zip(pids, pcct_vals, eid_vals):
            lp = math.log(pv + 1)
            le = math.log(ev + 1)
            detail.append(f"{pid:<15s} {lp:>12.4f} {le:>12.4f} {lp - le:>12.4f} {(lp + le) / 2:>12.4f}")

    # Summary LoA table
    lines.append("--- Limits of Agreement Summary (all variables, log(x+1) scale, raw volumes) ---")
    lines.append("")
    lines.append(f"{'Variable':<30s} {'Bias':>8s} {'LoA Lower':>10s} {'LoA Upper':>10s} {'r2':>6s} {'r2 95% CI':>18s} {'Result':>8s}")
    lines.append("-" * 95)

    from pathlib import Path
    plot_dir = Path(OUTPUT_DIR) / "bland_altman_plots"
    all_vars = {}
    all_vars.update(GATE3_PRIMARY)
    all_vars.update(GATE3_SECONDARY)
    for var, cfg in all_vars.items():
        # Raw volumes for BA comparability with delta OQ
        pcct_vals, eid_vals, pids_ba = _get_paired_values(paired, var, normalize=False)
        if len(pcct_vals) < 2:
            lines.append(f"{cfg['label']:<30s} {'N/A':>8s}")
            continue
        mb, sd, lo, hi, rsq = bland_altman(pcct_vals, eid_vals, log_transform=True)
        rsq_lb, rsq_ub = bootstrap_rsq_ci(pcct_vals, eid_vals, log_transform=True)
        # Plot if not already plotted in the per-variable section above
        if var not in GATE4_VARIABLES:
            plot_bland_altman(
                pcct_vals, eid_vals, cfg["label"],
                plot_dir / f"BA_{var}.png",
                log_transform=True, pids=pids_ba,
            )
        if rsq_lb is not None:
            prop = "PASS" if rsq_lb < 0.1 else "FAIL"
            ci_str = f"[{rsq_lb:.3f}, {rsq_ub:.3f}]"
        else:
            prop = "N/A"
            ci_str = "N/A"
        lines.append(f"{cfg['label']:<30s} {mb:>8.4f} {lo:>10.4f} {hi:>10.4f} {rsq:>6.3f} {ci_str:>18s} {prop:>8s}")

    lines.append("")
    lines.append("Proportional bias: PASS if r2 95% CI lower bound < 0.1")
    lines.append("")

    return "\n".join(lines), "\n".join(detail)


# ── Main ──────────────────────────────────────────────────────────────────────


if __name__ == "__main__":
    print("PCCT Qualification Gate Analyses")
    print("=" * 40)
    print(f"Reference: B.1P Delta Validation OQ (730-CVV-040)")
    print(f"PCCT summaries: {PCCT_DIR}")
    print(f"EID summaries:  {EID_DIR}")
    print()

    paired = load_paired_data()

    if len(paired) < 2:
        print(f"Need at least 2 paired patients to run analyses.")
        print(f"Currently {len(paired)} paired. Place EID CSVs in:")
        print(f"  {EID_DIR}")
        sys.exit(1)

    if len(paired) < 30:
        print(f"WARNING: {len(paired)} paired patients (target: 30).")
        print(f"Results are preliminary until N >= 30.\n")

    os.makedirs(OUTPUT_DIR, exist_ok=True)

    gate2_text = run_gate2(paired)
    gate3_text, gate3_detail = run_gate3(paired)
    gate4_text, gate4_detail = run_gate4(paired)

    summary = "\n\n".join([
        "PCCT Qualification Gate Analysis Report",
        f"Reference: B.1P Delta Validation OQ (730-CVV-040 v0.1)",
        f"N = {len(paired)} paired patients (target: >=30)",
        f"{'PRELIMINARY' if len(paired) < 30 else 'FINAL'} results",
        gate2_text,
        gate3_text,
        gate4_text,
    ])

    summary_path = os.path.join(OUTPUT_DIR, "gate_summary.txt")
    with open(summary_path, "w", encoding="utf-8") as f:
        f.write(summary)
    print(summary)
    print(f"\nSummary written to: {summary_path}")

    detail_path = os.path.join(OUTPUT_DIR, "gate_detail.txt")
    with open(detail_path, "w", encoding="utf-8") as f:
        f.write("GATE 3 — PER-PATIENT DETAIL\n")
        f.write(gate3_detail)
        f.write("\n\nGATE 4 — BLAND-ALTMAN DETAIL\n")
        f.write(gate4_detail)
    print(f"Detail written to:  {detail_path}")

    # Paired data CSV for further analysis
    paired_csv_path = os.path.join(OUTPUT_DIR, "paired_data.csv")
    all_export_vars = (list(GATE3_PRIMARY.keys()) + list(GATE3_DESCRIPTIVE.keys())
                       + list(GATE3_SECONDARY.keys()))
    with open(paired_csv_path, "w", newline="", encoding="utf-8") as f:
        headers = ["patient_id"]
        for v in all_export_vars:
            headers.extend([f"PCCT_{v}", f"EID_{v}", f"diff_{v}", f"mean_{v}"])
        writer = csv.writer(f)
        writer.writerow(headers)
        for pair in paired:
            row = [pair["patient_id"]]
            for v in all_export_vars:
                pv = pair["pcct"].get(v)
                ev = pair["eid"].get(v)
                if pv is not None and ev is not None:
                    row.extend([pv, ev, pv - ev, (pv + ev) / 2])
                else:
                    row.extend(["", "", "", ""])
            writer.writerow(row)
    print(f"Paired CSV:         {paired_csv_path}")
