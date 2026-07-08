# Gate analysis variants (N=25, latest 2026-07-07 workitem data)

Each file is the full `gate_summary.txt` from `run_gate_analyses.py` under a
different estimator/criterion configuration, so the choices can be compared
side by side. All use the same paired data; `gate_detail.txt`, `paired_data.csv`,
and the Bland-Altman figures are config-invariant (they don't depend on the
wCV/bias switches). The repo's top-level `gate_results/gate_summary.txt` matches
**variant B** (the plain default).

| File | Command | wCV estimator | scanner term | Gate 4 bias criterion |
|---|---|---|---|---|
| `summary_A_legacy-rmsrel.txt` | `--wcv-method rms-rel` | legacy rms-rel | no | pct-threshold |
| `summary_B_variance-component_DEFAULT.txt` | *(defaults)* | variance-component (Quan-Shih / OQ) | no | pct-threshold |
| `summary_C_variance-component_scanner-term.txt` | `--scanner-term` | variance-component | yes | pct-threshold |
| `summary_D_oq-bias-criterion.txt` | `--bias-criterion oq-ci-overlap` | variance-component | no | oq-ci-overlap |

See `../../tracker/statistical-methodology.md` for the rationale, the
verdict-impact tables, and caveats. Key points:
- Legacy → variance-component corrects the wCV (esp. the log branch); flips
  Lumen/Vessel process outputs PASS→FAIL on the canonical region.
- Scanner term removes the systematic modality bias from the wCV → process
  outputs come within OQ (marginal).
- oq-ci-overlap replaces the project-specific `|bias|<X%` Gate 4 threshold with
  95% CI overlap vs 730-CVV-040 Table 6 (OQ-consistent). NonCALC Matrix and Total
  Plaque bias FAIL under both criteria (real modality bias).

**Caveats:** N=25 preliminary (target ≥30); the SUB-SEGMENT section in every file
is on stale (pre-07-07) segmentations and needs regeneration; process-output bias
has no OQ reference table (descriptive under oq-ci-overlap).
