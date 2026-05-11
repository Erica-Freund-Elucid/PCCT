"""Generate Bland-Altman plots for the sensitivity analysis (excl PT-119, PT-133, PT-142).

Outputs to gate_results/bland_altman_plots/sensitivity_excl_119_133_142/.
Reads paired_data.csv (vessel-overlap totals from canonical run_gate_analyses.py run).
"""
import csv
import os
import sys

# Use the plot helpers from run_gate_analyses.py
PCCT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))
sys.path.insert(0, PCCT_ROOT)
from run_gate_analyses import plot_bland_altman, GATE4_VARIABLES, GATE3_PRIMARY, GATE3_SECONDARY

EXCLUDE = {"PT-119", "PT-133", "PT-142"}
PAIRED_CSV = os.path.join(PCCT_ROOT, "gate_results", "paired_data.csv")
OUT_DIR = os.path.join(
    PCCT_ROOT, "gate_results", "bland_altman_plots",
    "sensitivity_excl_119_133_142",
)
os.makedirs(OUT_DIR, exist_ok=True)


def main():
    rows = list(csv.DictReader(open(PAIRED_CSV, encoding="utf-8")))
    rows = [r for r in rows if r["patient_id"] not in EXCLUDE]
    print(f"N after exclusion: {len(rows)}")

    process_vars = set(GATE3_PRIMARY.keys())
    plaque_vars = set(GATE3_SECONDARY.keys())

    all_vars = {}
    all_vars.update(GATE3_PRIMARY)
    all_vars.update(GATE3_SECONDARY)

    for var, cfg in all_vars.items():
        pcct_col = f"PCCT_{var}"
        eid_col = f"EID_{var}"
        pids, pccts, eids = [], [], []
        for r in rows:
            try:
                p = float(r[pcct_col]); e = float(r[eid_col])
                pids.append(r["patient_id"]); pccts.append(p); eids.append(e)
            except (ValueError, KeyError):
                pass
        if len(pccts) < 2:
            print(f"  {var}: insufficient data ({len(pccts)})"); continue

        log_t = (var in process_vars)
        ref_b = ref_l = None
        if var in GATE4_VARIABLES and log_t:
            ref_b = GATE4_VARIABLES[var].get("ref_bias")
            ref_l = GATE4_VARIABLES[var].get("ref_loa")

        out_path = os.path.join(OUT_DIR, f"BA_{var}.png")
        from pathlib import Path
        plot_bland_altman(
            pccts, eids, cfg["label"],
            Path(out_path),
            log_transform=log_t,
            ref_bias=ref_b, ref_loa=ref_l,
            pids=pids,
        )
        print(f"  wrote {os.path.basename(out_path)} (scale={'log' if log_t else 'un'}, N={len(pccts)})")


if __name__ == "__main__":
    main()
