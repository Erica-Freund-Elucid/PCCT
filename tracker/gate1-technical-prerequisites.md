# Gate 1 — Technical & Image Quality Prerequisites

> **Required** — Scanner not eligible without passing all items

---

### 1.1 DICOM Compliance

- **Threshold:** Full DICOM 3.0; mandatory tags present
- **Method:** Automated tag validation on ingestion
- **Owner:**
- **Status:** Not started
- **Evidence:**
- **Notes:**

---

### 1.2 Contrast Timing

- **Threshold:** Peak aortic HU ≥ 300 in ≥ 90% of cases
- **Method:** Aortic ROI measurement across paired dataset
- **Owner:**
- **Status:** Not started
- **Evidence:**
- **Notes:**

---

### 1.3 Image Noise (SNR/CNR)

- **Threshold:** Non-inferior to reference CTA (±15%)
- **Method:** Uniform phantom or aortic ROI SD comparison
- **Owner:**
- **Status:** Not started
- **Evidence:**
- **Notes:**

---

### 1.4 Reconstruction Kernel

- **Threshold:** Soft-tissue kernel available; sharp kernel documented
- **Method:** Protocol review + DICOM header
- **Owner:** Murali
- **Status:** Complete
- **Evidence:** https://elucidbio.sharepoint.com/:w:/r/Technology/Documents/IP/PCCT/PT142_combined_analysis_COMPLETE_FINAL_v1.0.docx?d=w0fe5c7e3207c4790930cac6308e7bde7&csf=1&web=1&e=mQv35o
- **Notes:** Kernel family influences local vessel segment and patient level plaque characterization. Recon settings shoud be standardized for subsequent evaluations especially reproducibility and bias-agreement analyses
