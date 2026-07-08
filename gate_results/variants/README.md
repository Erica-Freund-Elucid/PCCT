# Gate analysis variants (N=25, latest 2026-07-07 workitem data)

Each file is the full `gate_summary.txt` from `run_gate_analyses.py` under a
different estimator/criterion configuration, so the choices can be compared
side by side. All use the same paired data; `gate_detail.txt`, `paired_data.csv`,
and the Gate-4 Bland-Altman figures are config-invariant (they don't depend on
the wCV switches). The repo's top-level `gate_results/gate_summary.txt` matches
**variant B** (the default).

| File | Command | wCV estimator | scanner term | Gate 4 bias criterion |
|---|---|---|---|---|
| `summary_A_legacy-rmsrel.txt` | `--wcv-method rms-rel --no-scanner-term` | legacy rms-rel | no | pct-threshold |
| `summary_B_default_var-comp_scanner-term.txt` | *(defaults)* | variance-component (Quan-Shih / OQ) | **yes (default)** | pct-threshold |
| `summary_C_no-scanner-term.txt` | `--no-scanner-term` | variance-component | no | pct-threshold |
| `summary_D_oq-bias-criterion.txt` | `--bias-criterion oq-ci-overlap` | variance-component | yes | oq-ci-overlap |

See `../../tracker/statistical-methodology.md` for the rationale, the
verdict-impact tables, and caveats. Key points:
- **Default (B)** is variance-component wCV **with the scanner term** — the
  like-for-like basis for the Gate 3 wCV-vs-OQ acceptance (OQ limit is a bias-free
  inter-operator dispersion, so the systematic scanner bias is excluded here and
  assessed in Gate 4). Under it all 7 endpoints' PCCT wCV CIs overlap the delta OQ.
- **B vs C** shows the scanner-term effect: without it (C) the raw cross-scanner
  wCV is inflated by the systematic modality bias and Lumen/Vessel/Wall fail the
  canonical CI-overlap.
- **A** (legacy rms-rel) is the original method/behavior, retained for traceability.
- **D** uses the OQ-consistent Gate 4 bias criterion: NonCALC Matrix and Total
  Plaque bias FAIL (real modality bias).

**Caveats:** N=25 preliminary (target ≥30); the SUB-SEGMENT section in every file
is on stale (pre-07-07) segmentations and needs regeneration; process-output bias
has no OQ reference table (descriptive under oq-ci-overlap).
