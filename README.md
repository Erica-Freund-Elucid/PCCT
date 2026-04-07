# PCCT Scanner Qualification Framework

Qualification criteria for Photon-Counting CT (PCCT) scanner validation against reference CTA. A new scanner must pass all four gates for final sign-off.

| Label | Meaning |
|---|---|
| **Required** | Gate criterion — must pass |
| **Performance** | Quantitative benchmark |
| **Advisory** | Recommended review |
| **Workflow step** | Process check |

---

## 1 · Technical & Image Quality Prerequisites

> **Required** — Gate 1

| Criterion | Threshold | Method |
|---|---|---|
| DICOM compliance | Full DICOM 3.0; mandatory tags present | Automated tag validation on ingestion |
| Contrast timing | Peak aortic HU ≥ 300 in ≥ 90% of cases | Aortic ROI measurement across paired dataset |
| Image noise (SNR/CNR) | Non-inferior to reference CTA (±15%) | Uniform phantom or aortic ROI SD comparison |
| Reconstruction kernel | Soft-tissue kernel available; sharp kernel documented | Protocol review + DICOM header |

---

## 2 · Workflow Integration Checks

> **Workflow step** — Gate 2

| Step | Pass criterion | Notes |
|---|---|---|
| DICOM ingestion | 100% successful ingestion; zero WIID assignment failures | Run full paired dataset through ingest pipeline |
| Centerline / vessel tree | Extraction success rate ≥ 95%; major vessel coverage equivalent to reference | Compare vessel segment count and length distribution |
| Lumen & wall init | Auto-initialization without manual override in ≥ 85% of segments | Log override frequency per scanner type |
| Lumen & wall editing | Edit rate not significantly higher than reference CTA (use 95% CI overlap as threshold) | Track editing events per case |
| Plaque quantification | See section 3 — core quantitative criteria | |
| Report generation | All required fields populated; no missing-data warnings | ECR Review in physician viewer and pdf reports for N ≥ 10 cases |

---

## 3 · Core Quantitative Reproducibility — Paired CTA vs. PCCT

> **Performance** — Gate 3

Assessed at patient level using wCV (within-subject coefficient of variation) of un-edited output from paired CTA/PCCT acquisitions of the same patient. Manual centerline definition for N ≥ 30 cases preferred to isolate variability due to image modality alone. Allows direct comparison to validated OQ. New scanner must meet **all** thresholds to pass.

| Measure | Reference wCV (existing) | New scanner threshold | Rationale |
|---|---|---|---|
| Lumen volume and vessel Length | < 10% wCV | ≤ 10% wCV (non-inferiority) | Tightest tolerance — primary anatomical endpoint |
| Calcified plaque volume | < 20% wCV | ≤ 20% wCV | High-attenuation segmentation; threshold-sensitive |
| Wall volume | < 30% wCV | ≤ 30% wCV | Outer wall boundary more variable; softer tolerance |
| Total plaque volume | < 30% wCV | ≤ 30% wCV | Composite of all plaque components |

### Statistical design

- Minimum **N = 30** paired patients (same-day or same-week CTA + PCCT)
- Compute wCV = (SD of differences / mean) × 100 per patient, then average across cohort
- Report **95% CI** on wCV
- If CI upper bound exceeds threshold, the criterion **fails** even if point estimate passes

---

## 4 · Systematic Bias & Agreement

> **Performance** — Gate 4

| Measure | Threshold | Method |
|---|---|---|
| Mean bias — lumen volume | \|bias\| < 5% of mean lumen volume | Bland-Altman plot; report LoA |
| Mean bias — calcified plaque | \|bias\| < 10% of mean calc volume | Bland-Altman plot; proportional bias test |
| Limits of agreement | LoA within ±1.96 SD; no proportional bias (r² < 0.1) | Regress residuals on mean |

---

## 5 · Contextual & Operational Checks

> **Advisory**

| Check | Recommendation |
|---|---|
| HU calibration stability | Verify calcified plaque HU threshold (typically 130 HU) remains valid on new scanner. If PCCT uses spectral data, map to conventional HU equivalent before thresholding. |
| Radiation dose | Document DLP/CTDI; confirm scanner does not require dose increases that offset clinical utility |
| Patient subgroup stratification | Report wCV separately for high vs. low calcium burden and BMI to check for heterogeneity |
| Reader/operator variability | Confirm inter-reader editing variability is within expected range on new scanner images |
| Compatability with ongoing projects | Confirm no _new_ failure modes with AVTE/AWAL |
| Software version lock | Qualify against software version where plaque algos equivalent to Nov2025 or later; re-qualification required for major algorithm updates |

---

## Qualification Decision Gate

| Gate | Requirement | Outcome on failure |
|---|---|---|
| **Gate 1** — Technical prerequisites | All section 1 items pass | Scanner not eligible; re-test with protocol adjustments |
| **Gate 2** — Workflow integration | All section 2 items pass | Engineering remediation required before performance testing |
| **Gate 3** — Reproducibility | All four wCV thresholds met (section 3) | Scanner conditionally failed; may expand N or adjust acquisition protocol and re-test once |
| **Gate 4** — Bias & agreement | All section 4 items pass | Flag for clinical review; may restrict use to specific indications |

**Final sign-off requires passing all four gates.** Advisory checks (section 5) inform labelling and use restrictions but do not block qualification.
