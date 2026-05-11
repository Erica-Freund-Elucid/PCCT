# Gate 2 -- Workflow Integration Checks

> **Workflow step** -- Engineering remediation required before performance testing

---

### 2.1 DICOM Ingestion

- **Threshold:** 100% successful ingestion; zero WIID assignment failures
- **Method:** Run full paired dataset through ingest pipeline
- **Owner:**
- **Status:** In progress
- **Evidence:**
- **Notes:** 30 PCCT and 30 EID workitems being tracked per workflow spreadsheet (excl. PT-149 EID rejected, PT-120 EID NA). 28 PCCT and 26 EID have summary CSVs. PT-124 paired but excluded from Gate 3/4: PCCT processed LeftCoronary, EID processed RightCoronary -- no overlap. Partial-overlap pairs analyzed on intersection only: PT-131, PT-136, PT-142, PT-150, PT-158, PT-161. PT-152 PCCT (wi-bb33e995) pending -- summary CSV not yet generated. PT-163 held for 1024x1024 separate evaluation. PT-152 EID, PT-158 EID Aorta.nrrd missing -- excluded from SNR.

---

### 2.2 Centerline / Vessel Tree

- **Threshold:** Extraction success rate >= 95%; major vessel coverage equivalent to reference
- **Method:** Compare vessel segment count and length distribution between PCCT and EID
- **Owner:**
- **Status:** REVIEW
- **Evidence:** gate_results/gate_summary.txt (Gate 2 section)
- **Notes:** N=25 paired patients (vessel-overlap). Mean length difference +26.3% on overlapping vessels (down from +44% with target-overlap-all-vessels). Confirms most of the apparent length disparity came from PCCT tracing extra distal vessels (PDA, PLB, marginals, diagonals) that EID did not. On the shared anatomy, length agreement is reasonable.

---

### 2.3 Lumen & Wall Initialization

- **Threshold:** Auto-initialization without manual override in >= 85% of segments
- **Method:** Log override frequency per scanner type
- **Owner:**
- **Status:** Not started
- **Evidence:**
- **Notes:**

---

### 2.4 Lumen & Wall Editing

- **Threshold:** Edit rate / analyst effort not significantly higher than reference CTA
- **Method:** Per-case analyst effort scored on 1-5 scale (1=most effort, 5=least). Tracked alongside lumen/wall/tissue-composition quality and app-stability incidents.
- **Owner:**
- **Status:** PASS (preliminary; PCCT and CCTA effort effectively equivalent on paired cases)
- **Evidence:** PCCT_CCTA_Case_Summaries.xlsx (Case Summaries sheet)
- **Notes:**
  - Per-scan summary: PCCT mean effort 2.68 (median 3.0) across N=31; CCTA mean effort 2.65 (median 3.0) across N=26.
  - Paired comparison restricted to the Gate 3/4 cohort (N=23 patients with numeric effort on both scans -- PT-156 and PT-158 CCTA effort "not recorded" in spreadsheet):
    - PCCT mean 2.85, CCTA mean 2.74; delta (PCCT - CCTA) mean +0.11, median 0
    - PCCT easier than CCTA in 7/23 (30%); PCCT harder in 9/23 (39%); same in 7/23 (30%)
  - Lumen Quality (Good/Fair/Bad): PCCT 9/9/13 vs CCTA 5/12/12 -- PCCT marginally more "Good" lumens.
  - Wall Quality: PCCT 9/11/11 vs CCTA 5/13/10 -- comparable.
  - App Stability: both scan types had Save & Generate Analysis crash incidents -- engineering issue, not a modality effect.
  - **Conclusion:** PCCT does not impose meaningfully higher analyst effort than the reference CTA workflow. Mean delta (+0.11) is well within the 0.5-point score granularity; patient-level split is roughly symmetric (7 easier / 9 harder / 7 same).

---

### 2.5 Plaque Quantification

- **Threshold:** See Gate 3 -- core quantitative criteria
- **Method:** Paired CTA vs. PCCT comparison
- **Owner:**
- **Status:** In progress (N=25 of 30 paired, vessel-overlap)
- **Evidence:** gate_results/gate_summary.txt (Gate 3 section)
- **Notes:** **All 7 endpoints (Lumen, Wall, Vessel, CALC, LRNC, NonCALC, Total Plaque) PASS** under vessel-overlap normalization (canonical as of 2026-05-08). Wall recovered from FAIL by switching from target-overlap to vessel-overlap.

---

### 2.6 Report Generation

- **Threshold:** All required fields populated; no missing-data warnings
- **Method:** ECR Review in physician viewer and pdf reports for N >= 10 cases
- **Owner:**
- **Status:** Not started
- **Evidence:**
- **Notes:**
