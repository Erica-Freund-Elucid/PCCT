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
| **Gate 1** -- Technical prerequisites | All items pass | Scanner not eligible; re-test with protocol adjustments |
| **Gate 2** -- Workflow integration | All items pass | Engineering remediation required before performance testing |
| **Gate 3** -- Reproducibility | 95% CI overlap with delta OQ | Scanner conditionally failed; may expand N or adjust protocol |
| **Gate 4** -- Bias & agreement | Bland-Altman within reference LoA | Flag for clinical review; may restrict use to specific indications |

**Final sign-off requires passing all four gates.** Advisory checks inform labelling and use restrictions but do not block qualification.

---

## 1 - Technical & Image Quality Prerequisites (Gate 1)

| Criterion | Threshold | Method | Status |
|---|---|---|---|
| DICOM compliance | Full DICOM 3.0; mandatory tags present | Automated tag validation on ingestion | Complete |
| Contrast timing | Peak aortic HU >= 300 in >= 90% of cases | Aortic ROI (10mm cube at centroid) | PCCT PASS (93%), EID FAIL (44%) |
| Image noise (SNR) | Non-inferior to EID (ratio >= 0.85) | Aortic ROI SNR = mean_HU / std_HU | PASS (ratio 2.01, PCCT 2x better) |
| Reconstruction kernel | Soft-tissue kernel available | Protocol review + DICOM header | Complete |

PCCT: mean SNR 18.5, mean noise 21.5 HU (N=28). EID: mean SNR 9.7, mean noise 32.8 HU (N=9). PCCT lower noise in 8/9 paired cases.

---

## 2 - Workflow Integration Checks (Gate 2)

| Step | Pass criterion | Status |
|---|---|---|
| DICOM ingestion | 100% success; zero WIID failures | In progress (28/30 PCCT, 9/13 EID) |
| Centerline / vessel tree | >= 95% extraction; coverage equivalent | REVIEW (PCCT traces longer vessels) |
| Lumen & wall initialization | Auto-init >= 85% of segments | Not started |
| Lumen & wall editing | Edit rate not significantly higher than reference | Not started |
| Plaque quantification | See Gate 3 | In progress |
| Report generation | All fields populated; no missing-data warnings | Not started |

---

## 3 - Core Quantitative Reproducibility (Gate 3)

**Acceptance:** 95% CI overlap between PCCT wCV and delta OQ wCV (730-CVV-040).

**Method:** Volumes length-normalized, then log(x+1) transform. wCV = sqrt(mean(d^2/(2m^2))) x 100. Bootstrap 95% CIs. N >= 30 paired patients required.

### Primary Endpoints

| Variable | PCCT log-wCV [95% CI] | Delta OQ [95% CI] | Original OQ [95% CI] | Status |
|---|---|---|---|---|
| Lumen Volume | 12.97% [8.98%, 16.21%] | 7.5% [5.9%, 9.3%] | 7.64% [5.52%, 8.55%] | **PASS** |
| Wall Volume | 24.91% [13.78%, 35.14%] | 13.1% [10.2%, 16.0%] | 9.51% [6.60%, 11.16%] | **PASS** |
| Vessel Volume | 12.44% [8.48%, 15.73%] | 8.7% [6.8%, 10.6%] | 6.03% [4.44%, 6.74%] | **PASS** (log) |

### Secondary Endpoints (plaque volumes)

| Variable | PCCT log-wCV [95% CI] | Delta OQ [95% CI] | Original OQ | Status |
|---|---|---|---|---|
| CALC Volume | 26.17% [13.27%, 35.68%] | 4.5% [3.5%, 5.4%] | 13.9% | FAIL |
| LRNC Volume | 80.52% [54.67%, 104.07%] | 5.4% [4.1%, 6.9%] | 59.3% | FAIL |
| NonCALC Matrix | 35.22% [19.69%, 48.21%] | 13.6% [10.8%, 16.4%] | 58.6% | FAIL |
| Total Plaque | 27.37% [16.16%, 38.63%] | 13.2% [10.4%, 16.0%] | 44.1% | FAIL |

### Descriptive

| Variable | PCCT untransformed wCV | Delta OQ | Note |
|---|---|---|---|
| Vessel Length | 44.37% [28.94%, 55.46%] | 13.88% [10.72%, 17.38%] | Not acceptance criterion (730-CVV-040) |

*All results preliminary at N=9. CIs will narrow substantially at N=30.*

---

## 4 - Systematic Bias & Agreement (Gate 4)

**Method:** Bland-Altman on log(x+1) scale, raw volumes (not length-normalized) for comparability with delta OQ (730-CVV-040 Table 6). Proportional bias r^2 with bootstrap 95% CI.

| Variable | Log bias | LoA | r^2 [95% CI] | Bias result | Prop. bias |
|---|---|---|---|---|---|
| Lumen | 0.24 | [-0.29, 0.77] | 0.013 [0.000, 0.661] | FAIL (21%) | **PASS** |
| Wall | 0.21 | [-0.51, 0.93] | 0.140 [0.003, 0.672] | -- | **PASS** |
| Vessel | 0.21 | [-0.17, 0.60] | 0.117 [0.001, 0.605] | -- | **PASS** |
| CALC | 0.52 | [-0.35, 1.38] | 0.090 [0.003, 0.567] | FAIL (43%) | **PASS** |
| LRNC | 0.86 | [-1.45, 3.17] | 0.055 [0.000, 0.855] | FAIL (65%) | **PASS** |
| NonCALC Matrix | 0.08 | [-0.80, 0.97] | 0.225 [0.019, 0.808] | **PASS** (1%) | **PASS** |
| Total Plaque | 0.12 | [-0.38, 0.62] | 0.188 [0.005, 0.765] | FAIL (18%) | **PASS** |

All 7 variables pass the proportional bias test (r^2 CI lower bound < 0.1). Lumen bias driven by PCCT having longer vessel tracing extent, not segmentation difference.

---

## 5 - Advisory Checks

| Check | Recommendation | Status |
|---|---|---|
| HU calibration stability | Verify 130 HU CALC threshold valid on PCCT | Not started |
| Radiation dose | Document DLP/CTDI | Not started |
| Subgroup stratification | Report wCV by calcium burden and BMI | Not started |
| Reader/operator variability | Inter-reader editing variability on PCCT images | Not started |
| Compatibility with ongoing projects | No new failure modes with AVTE/AWAL | Not started |
| Software version lock | Qualify against Nov2025+ plaque algorithms | Not started |

---

## Scripts

| Script | Description |
|---|---|
| `run_gate_analyses.py` | Main analysis: Gates 1-4 from paired PCCT/EID summary CSVs + SNR data |
| `generate_tracker.py` | Generates Excel tracker from markdown gate files |

### Usage

```
# Place PCCT summary CSVs in workitem_summaries/PCCT/
# Place EID summary CSVs in workitem_summaries/EID/
# Place SNR results in gate_results/snr_pcct.csv and snr_eid.csv
python run_gate_analyses.py
```

### Outputs (gate_results/)

| File | Description |
|---|---|
| `gate_summary.txt` | Full report with pass/fail for each criterion |
| `gate_detail.txt` | Per-patient breakdowns for each variable |
| `paired_data.csv` | Tidy paired data for further analysis |
| `bland_altman_plots/` | Bland-Altman plots per variable |
| `snr_pcct.csv` / `snr_eid.csv` | Gate 1 SNR measurements from aortic ROI |

---

## Tracker

Each criterion has a dedicated entry in `tracker/` with **Owner**, **Status**, **Evidence**, and **Notes**.

| File | Contents |
|---|---|
| [Gate 1](tracker/gate1-technical-prerequisites.md) | DICOM compliance, contrast timing, noise, kernel |
| [Gate 2](tracker/gate2-workflow-integration.md) | Ingestion, centerline, lumen/wall, reporting |
| [Gate 3](tracker/gate3-reproducibility.md) | wCV for lumen, wall, vessel vol, plaque volumes |
| [Gate 4](tracker/gate4-bias-agreement.md) | Bland-Altman, LoA, proportional bias |
| [Advisory](tracker/advisory-operational-checks.md) | HU calibration, dose, subgroups, reader variability |
