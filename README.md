# PCCT Scanner Qualification Framework

Qualification criteria for Photon-Counting CT (PCCT) scanner validation against reference EID CTA. Acceptance is based on non-inferiority to the validated B.1P inter-operator reproducibility established in the original OQ (4-B1P-033) and delta validation OQ (730-CVV-040).

## Reference Documents

| Document | Description |
|---|---|
| 4-B1P-033 v2.0 | B.1P Operational Qualification Report (original clearance) |
| 730-CVV-040 v0.1 | B.1P Delta Validation OQ Report |
| 730-CVV-017 | B.1P November 2024 Release NC Matrix Validation Report |

## Gate Structure

| Gate | Requirement | Outcome on failure |
|---|---|---|
| **Gate 1** -- Technical prerequisites | All section 1 items pass | Scanner not eligible; re-test with protocol adjustments |
| **Gate 2** -- Workflow integration | All section 2 items pass | Engineering remediation required before performance testing |
| **Gate 3** -- Reproducibility | 95% CI overlap with delta OQ | Scanner conditionally failed; may expand N or adjust acquisition protocol |
| **Gate 4** -- Bias & agreement | Bland-Altman within reference LoA | Flag for clinical review; may restrict use to specific indications |

**Final sign-off requires passing all four gates.** Advisory checks inform labelling and use restrictions but do not block qualification.

---

## 3 - Core Quantitative Reproducibility (Gate 3)

**Acceptance criterion:** 95% CI overlap between PCCT qualification wCV and B.1P Delta Validation OQ wCV (730-CVV-040). Cross-scanner variability must be non-inferior to validated inter-operator variability.

**Primary metric:** Log-wCV for volume endpoints (measurement error proportional to magnitude). Untransformed wCV for vessel length (additive error, descriptive only per 730-CVV-040).

**Volumes are length-normalized** before computing wCV to isolate variability due to image modality from vessel extent differences.

### Primary Endpoints

| Variable | Delta OQ log-wCV (95% CI) | Original OQ log-wCV (95% CI) |
|---|---|---|
| Lumen Volume | 7.5% [5.9%, 9.3%] | 7.64% [5.52%, 8.55%] |
| Wall Volume | 13.1% [10.2%, 16.0%] | 9.51% [6.60%, 11.16%] |
| Vessel Volume | 8.7% [6.8%, 10.6%] | 6.03% [4.44%, 6.74%] |

### Secondary Endpoints (plaque volumes)

| Variable | Delta OQ log-wCV (95% CI) | Original OQ (untransformed) |
|---|---|---|
| CALC Volume | 4.5% [3.5%, 5.4%] | 13.9% |
| LRNC Volume | 5.4% [4.1%, 6.9%] | 59.3% |
| NonCALC Matrix Volume | 13.6% [10.8%, 16.4%] | 58.6% (Total NCP) |
| Total Plaque Volume | 13.2% [10.4%, 16.0%] | 44.1% |

### Descriptive (not acceptance criterion)

| Variable | Delta OQ untransformed wCV (95% CI) | Original OQ log-wCV (95% CI) |
|---|---|---|
| Vessel Length | 13.88% [10.72%, 17.38%] | 5.27% [3.43%, 6.52%] |

Vessel length is descriptive only per 730-CVV-040 recommendation. Length variability reflects analyst termination judgment, not segmentation error.

### Statistical Design

- Minimum **N = 30** paired patients (same-patient PCCT + EID CTA)
- wCV = sqrt(mean(d^2 / (2m^2))) x 100
- 95% CI via 2000-sample bootstrap
- log(x+1) transformation for volume endpoints
- Pass: PCCT 95% CI overlaps with delta OQ 95% CI

---

## 4 - Systematic Bias & Agreement (Gate 4)

**Method:** Bland-Altman on log(x+1) scale using raw (not length-normalized) volumes for direct comparability with delta OQ reference (730-CVV-040 Table 6).

| Measure | Threshold | Method |
|---|---|---|
| Mean bias -- lumen volume | \|bias\| < 5% of mean | Bland-Altman; report LoA |
| Mean bias -- plaque volumes | \|bias\| < 10% of mean | Bland-Altman; proportional bias test |
| Proportional bias | r^2 95% CI lower bound < 0.1 | Bootstrap CI on regression of residuals on mean |
| Limits of agreement | LoA within +/-1.96 SD | Compare with delta OQ reference LoA |

---

## Scripts

| Script | Description |
|---|---|
| `run_gate_analyses.py` | Main analysis: computes Gate 2, 3, 4 from paired PCCT/EID summary CSVs |
| `generate_tracker.py` | Generates Excel tracker from markdown gate files |

### Usage

```
# Place PCCT summary CSVs in workitem_summaries/PCCT/
# Place EID summary CSVs in workitem_summaries/EID/
python run_gate_analyses.py
```

### Outputs (gate_results/)

| File | Description |
|---|---|
| `gate_summary.txt` | Full report with pass/fail for each criterion |
| `gate_detail.txt` | Per-patient breakdowns for each variable |
| `paired_data.csv` | Tidy paired data for further analysis |
| `bland_altman_plots/` | Bland-Altman plots per variable |
| `snr_results.csv` | Gate 1.3 SNR measurements from aortic ROI |

---

## Tracker

Each criterion has a dedicated entry in the `tracker/` folder with fields for **Owner**, **Status**, **Evidence**, and **Notes**.

| File | Contents |
|---|---|
| [Gate 1 -- Technical Prerequisites](tracker/gate1-technical-prerequisites.md) | DICOM compliance, contrast timing, noise, kernel |
| [Gate 2 -- Workflow Integration](tracker/gate2-workflow-integration.md) | Ingestion, centerline, lumen/wall, reporting |
| [Gate 3 -- Reproducibility](tracker/gate3-reproducibility.md) | wCV for lumen, wall, vessel vol, plaque volumes |
| [Gate 4 -- Bias & Agreement](tracker/gate4-bias-agreement.md) | Bland-Altman, LoA, proportional bias |
| [Advisory Checks](tracker/advisory-operational-checks.md) | HU calibration, dose, subgroups, reader variability |
