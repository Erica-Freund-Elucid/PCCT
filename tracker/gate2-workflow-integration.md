# Gate 2 -- Workflow Integration Checks

> **Workflow step** -- Engineering remediation required before performance testing

---

### 2.1 DICOM Ingestion

- **Threshold:** 100% successful ingestion; zero WIID assignment failures
- **Method:** Run full paired dataset through ingest pipeline
- **Owner:**
- **Status:** In progress
- **Evidence:**
- **Notes:** 30 PCCT and 13 EID workitems assigned per workflow spreadsheet. 28 PCCT and 9 EID have summary CSVs (processing complete). 2 PCCT (PT-152, PT-163) and 4 EID (PT-124, PT-130, PT-142, PT-163) pending.

---

### 2.2 Centerline / Vessel Tree

- **Threshold:** Extraction success rate >= 95%; major vessel coverage equivalent to reference
- **Method:** Compare vessel segment count and length distribution between PCCT and EID
- **Owner:**
- **Status:** REVIEW
- **Evidence:** gate_results/gate_summary.txt (Gate 2 section)
- **Notes:** N=9 paired patients. Mean length difference +86.5% (PCCT longer than EID). Mean segment count difference +3.7 (PCCT more segments). Large differences in PT-129 (107%), PT-158 (177%), PT-165 (184%) suggest EID scans have significantly shorter vessel tracing -- likely due to lower contrast/image quality on Somatom Force. This is controlled by length normalization in Gate 3.

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

- **Threshold:** Edit rate not significantly higher than reference CTA (use 95% CI overlap as threshold)
- **Method:** Track editing events per case
- **Owner:**
- **Status:** Not started
- **Evidence:**
- **Notes:**

---

### 2.5 Plaque Quantification

- **Threshold:** See Gate 3 -- core quantitative criteria
- **Method:** Paired CTA vs. PCCT comparison
- **Owner:**
- **Status:** In progress (N=9 of 30 paired)
- **Evidence:** gate_results/gate_summary.txt (Gate 3 section)
- **Notes:** Primary endpoints (lumen, wall, vessel vol) passing CI overlap. Plaque volumes have wider CIs at N=9.

---

### 2.6 Report Generation

- **Threshold:** All required fields populated; no missing-data warnings
- **Method:** ECR Review in physician viewer and pdf reports for N >= 10 cases
- **Owner:**
- **Status:** Not started
- **Evidence:**
- **Notes:**
