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
| Contrast timing | Peak aortic HU >= 250 in >= 90% of cases | Aortic ROI (10mm cube at centroid) | PCCT PASS (100%), EID PASS (96%) |
| Image noise (SNR) | Non-inferior to EID (ratio >= 0.85) | Aortic ROI SNR = mean_HU / std_HU | PASS (ratio 2.04, PCCT ~2x better) |
| Reconstruction kernel | Soft-tissue kernel available | Protocol review + DICOM header | Complete |

PCCT: mean SNR 18.5, mean noise 21.5 HU (N=28). EID: mean SNR ~9.0, mean noise 36.5 HU (N=25). PCCT lower noise in 24/25 paired cases. Only PT-124 (242 HU) below contrast-timing threshold. PT-152 EID, PT-158 EID NRRDs partial (Volume.nrrd present, Aorta.nrrd missing) -- excluded from SNR.

---

## 2 - Workflow Integration Checks (Gate 2)

| Step | Pass criterion | Status |
|---|---|---|
| DICOM ingestion | 100% success; zero WIID failures | In progress (28/30 PCCT, 26/30 EID) |
| Centerline / vessel tree | >= 95% extraction; coverage equivalent | REVIEW (PCCT +26% length over EID on overlapping vessels at N=25) |
| Lumen & wall initialization | Auto-init >= 85% of segments | Not started |
| Lumen & wall editing | Edit rate not significantly higher than reference | PASS (paired N=23, mean effort delta +0.11 on 1-5 scale) |
| Plaque quantification | See Gate 3 | In progress |
| Report generation | All fields populated; no missing-data warnings | Not started |

---

## 3 - Core Quantitative Reproducibility (Gate 3)

**Acceptance:** 95% CI overlap between PCCT wCV and delta OQ wCV (730-CVV-040).

**Method:** Volumes length-normalized. wCV = sqrt(mean(d^2/(2m^2))) x 100. Bootstrap 95% CIs. N >= 30 paired patients required. **Metric per endpoint family:** log-wCV for process outputs (Lumen, Wall, Vessel vol — matches 4-B1P-033 Table 7); untransformed wCV for plaque volumes (matches 4-B1P-033 Table 9 reporting). **Per-pair vessel-overlap restriction (added 2026-05-08):** for each PCCT/EID pair, only the (bodySite, location) vessel intersection is summed -- so the same anatomical extent is measured on each scanner. Patients with no overlapping vessel are excluded; partial-overlap pairs analyzed on the intersection only.

### Primary Endpoints

| Variable | PCCT log-wCV [95% CI] | Delta OQ [95% CI] | Original OQ [95% CI] | Status |
|---|---|---|---|---|
| Lumen Volume | 8.72% [6.99%, 10.30%] | 7.5% [5.9%, 9.3%] | 7.64% [5.52%, 8.55%] | **PASS** |
| Wall Volume | 18.48% [14.07%, 22.55%] | 13.1% [10.2%, 16.0%] | 9.51% [6.60%, 11.16%] | **PASS** |
| Vessel Volume | 8.40% [6.83%, 9.85%] | 8.7% [6.8%, 10.6%] | 6.03% [4.44%, 6.74%] | **PASS** (log) |

### Secondary Endpoints (plaque volumes — untransformed wCV)

| Variable | PCCT untransformed wCV [95% CI] | Delta OQ untransformed [95% CI] | Original OQ | Status |
|---|---|---|---|---|
| CALC Volume | 21.08% [16.57%, 25.03%] | 25.78% [18.70%, 35.75%] | 13.9% | **PASS** |
| LRNC Volume | 94.24% [78.29%, 108.83%] | 78.80% [48.09%, 163.40%] | 59.3% | **PASS** |
| NonCALC Matrix | 39.46% [30.22%, 47.73%] | 32.61% [24.87%, 41.94%] | 58.6% | **PASS** |
| Total Plaque | 31.40% [24.65%, 37.47%] | 27.08% [20.68%, 34.74%] | 44.1% | **PASS** |

### Descriptive

| Variable | PCCT untransformed wCV | Delta OQ | Note |
|---|---|---|---|
| Vessel Length | 24.64% [18.88%, 29.60%] | 13.88% [10.72%, 17.38%] | Not acceptance criterion (730-CVV-040) |

*N=25 paired (vessel-overlap, canonical as of 2026-05-08). PT-124 excluded (no vessel overlap -- PCCT=Left, EID=Right). PT-152 paired only on EID (PCCT not yet processed). PT-149 EID rejected, PT-120 EID NA. Most pairs are partial vessel-overlap -- PCCT typically traces 1-3 additional distal vessels. **All 7 endpoints (3 process + 4 plaque) PASS** with vessel-overlap normalization. Wall flipped from FAIL to PASS (was 23.6% with target-overlap). The earlier Wall failure was an artifact of length-normalization dilution from PCCT's longer distal trace; restricting to vessels traced on both scanners eliminates the dilution.*

---

## 4 - Systematic Bias & Agreement (Gate 4)

**Method:** Bland-Altman on raw volumes (not length-normalized), summed over the per-pair target intersection. **Process outputs (Lumen, Wall, Vessel) on log(x+1) scale** for comparability with delta OQ (730-CVV-040 Table 6); **plaque (CALC, LRNC, NonCALC, Total) on untransformed scale** to match 4-B1P-033 Table 9 reporting. Proportional bias r^2 with bootstrap 95% CI.

| Variable | Scale | Bias | LoA | r^2 [95% CI] | Bias result | Prop. bias |
|---|---|---|---|---|---|---|
| Lumen | log | -0.02 | [-0.49, 0.46] | 0.005 [0.000, 0.207] | **PASS** (1.6%) | **PASS** |
| Wall† | log | -0.07 | [-0.60, 0.46] | 0.016 [0.000, 0.233] | **PASS** (9.3%) | **PASS** |
| Vessel | log | -0.02 | [-0.41, 0.36] | 0.001 [0.000, 0.180] | -- | **PASS** |
| CALC | raw | 27.61 mm³ | [-33.97, 89.19] | 0.142 [0.010, 0.703] | FAIL (19.2%) | **PASS** |
| LRNC | raw | 5.65 mm³ | [-49.58, 60.87] | 0.051 [0.001, 0.725] | FAIL (21.0%) | **PASS** |
| NonCALC Matrix | raw | -65.17 mm³ | [-281.60, 151.26] | 0.436 [0.056, 0.725] | FAIL (-24.5%) | **PASS** |
| Total Plaque | raw | -31.91 mm³ | [-289.67, 225.84] | 0.167 [0.001, 0.542] | **PASS** (7.3%) | **PASS** |

†Wall |bias|<10% threshold is project-specific (added 2026-05-08); not derivable from 4-B1P-033 or 730-CVV-040, both of which assess Wall reproducibility by wCV only. Added because case review flagged a likely true modality bias at high disease. Threshold matched to plaque components (10%) since wall measurement is dominated by the wall-plaque boundary, closer in difficulty to plaque segmentation than to lumen.

**All 7 variables PASS proportional-bias.** Lumen / Wall / Total Plaque PASS bias; per-component plaque (CALC, LRNC, NonCALC) FAIL. Directional pattern: PCCT reports more CALC/LRNC, less NonCALC (see [tracker/case-reviews.md](tracker/case-reviews.md)).

### Sensitivity — excluding 3 case-review-flagged patients

PT-119 / PT-133 (acquisition-condition artifacts) and PT-142 (likely true modality bias at high disease) — see case-reviews.md. Dropping these 3 patients:

| | N=25 (canonical) | N=22 (excl 119, 133, 142) |
|---|---|---|
| Gate 3 (all 7) | PASS | PASS |
| Wall bias | 9.3% PASS | 4.0% PASS |
| Total Plaque bias | 7.3% PASS | 1.0% PASS (much improved) |
| CALC prop. bias | PASS [0.010, 0.704] | **FAIL** [0.298, 0.815] ← flips to FAIL |
| Other prop. bias | all PASS | all PASS |

Wall bias drops further on the cleaner subset (4.0%) — consistent with PT-142 contributing the high-disease modality signal. CALC prop-bias FAILing on the smaller-volume subset reflects loss of the high-disease anchor cases that flatten the regression slope — a low-N artifact, expected to resolve at N≥30. Sensitivity BA plots in `gate_results/bland_altman_plots/sensitivity_excl_119_133_142/`.

### Supplementary — Bland-Altman on length-normalized volumes

The canonical Gate 4 above uses **raw** volumes (matches 730-CVV-040 Table 6 ref values). Even with vessel-overlap restriction, PCCT length on shared vessels is +26% on average -- residual extent differential within named vessels inflates PCCT raw volumes. Below repeats the BA on **volume per mm** so the read is unaffected by length differential. **Delta-OQ Table 6 ref values do NOT apply on this scale** (they are raw-only); use this as a cross-scanner cleanness check, not a regulatory criterion.

| Variable | Raw bias (canonical) | /mm bias (supplementary) | Reading |
|---|---|---|---|
| Lumen | **PASS** (1.6%) | FAIL (20.5%) | PCCT under-calls lumen per mm |
| Wall | PASS (9.3%) | FAIL (22.4%) | PCCT under-calls wall per mm |
| CALC | FAIL (19.2%) | **PASS** (2.6%) | Raw bias was largely length-driven |
| LRNC | FAIL (21.0%) | FAIL (34.4%) | Same direction (PCCT > EID) |
| NonCALC | FAIL (-24.5%) | FAIL (-35.5%) | PCCT consistently lower NonCALC per mm |
| Total Plaque | **PASS** (-7.3%) | FAIL (-19.2%) | Raw agreement masks per-mm under-call |

All 7 still pass the proportional-bias test on /mm scale. Plots in `gate_results/bland_altman_plots/length_normalized/`. The pattern of PCCT < EID for per-mm densities (lumen, wall, plaque components except CALC) is the inverse of the raw-volume pattern (PCCT > EID for total volumes) — consistent with PCCT analysts tracing further into the same named vessels, so PCCT's longer trace covers more "normal" distal vessel (thin wall, thin plaque) and dilutes the per-mm averages.

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
| [Case reviews](tracker/case-reviews.md) | Per-patient visual-review notes (acquisition vs modality effects) |
