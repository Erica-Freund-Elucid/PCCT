# Gate 2 — Workflow Integration Checks

> **Workflow step** — Engineering remediation required before performance testing

---

### 2.1 DICOM Ingestion

- **Threshold:** 100% successful ingestion; zero WIID assignment failures
- **Method:** Run full paired dataset through ingest pipeline
- **Owner:**
- **Status:** Not started
- **Evidence:**
- **Notes:**

---

### 2.2 Centerline / Vessel Tree

- **Threshold:** Extraction success rate ≥ 95%; major vessel coverage equivalent to reference
- **Method:** Compare vessel segment count and length distribution
- **Owner:**
- **Status:** Not started
- **Evidence:**
- **Notes:**

---

### 2.3 Lumen & Wall Initialization

- **Threshold:** Auto-initialization without manual override in ≥ 85% of segments
- **Method:** Log override frequency per scanner type
- **Owner:**
- **Status:** Not started
- **Evidence:**
- **Notes:**

---

### 2.4 Lumen & Wall Editing

- **Threshold:** Edit rate not significantly higher than reference CTA (use 95% CI overlap as threshold)
- **Method:** Track editing events per case, can we get this for ones Carolyn already did?
- **Owner:**
- **Status:** Not started
- **Evidence:**
- **Notes:**

---

### 2.5 Plaque Quantification

- **Threshold:** See Gate 3 — core quantitative criteria
- **Method:** Paired CTA vs. PCCT comparison
- **Owner:**
- **Status:** Not started
- **Evidence:**
- **Notes:**

---

### 2.6 Report Generation

- **Threshold:** All required fields populated; no missing-data warnings
- **Method:** ECR Review in physician viewer and pdf reports for N ≥ 10 cases use ones Carolyn already did
- **Owner:**
- **Status:** Not started
- **Evidence:**
- **Notes:**
