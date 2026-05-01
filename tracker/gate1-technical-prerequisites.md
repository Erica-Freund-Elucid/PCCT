# Gate 1 -- Technical & Image Quality Prerequisites

> **Required** -- Scanner not eligible without passing all items

---

### 1.1 DICOM Compliance

- **Threshold:** Full DICOM 3.0; mandatory tags present
- **Method:** Automated tag validation on ingestion
- **Owner:** Murali
- **Status:** Complete
- **Evidence:** All 30 PCCT cases (NAEOTOM Alpha) successfully ingested on ip3
- **Notes:**

---

### 1.2 Contrast Timing

- **Threshold:** Peak aortic HU >= 300 in >= 90% of cases
- **Method:** Aortic ROI measurement (10mm cube at aortic centroid via compute_snr.py)
- **Owner:**
- **Status:** PCCT PASS / EID FAIL
- **Evidence:** gate_results/snr_pcct.csv, gate_results/snr_eid.csv
- **Notes:** PCCT: 26/28 (93%) >= 300 HU -- PASS. EID: 4/9 (44%) >= 300 HU -- FAIL. EID failure driven by lower kVp protocols on Somatom Force (90-130 kVp). 5 EID cases below threshold: PT-140 (252 HU), PT-130 (266 HU), PT-116 (267 HU), PT-162 (279 HU), PT-129 (284 HU). This is a protocol difference, not a scanner limitation.

---

### 1.3 Image Noise (SNR/CNR)

- **Threshold:** Non-inferior to reference EID CTA (SNR ratio >= 0.85)
- **Method:** Aortic ROI SNR = mean_HU / std_HU (10mm cube at aortic centroid)
- **Owner:**
- **Status:** PASS
- **Evidence:** gate_results/snr_pcct.csv, gate_results/snr_eid.csv
- **Notes:** PCCT mean SNR 18.53 vs EID mean SNR 9.73. Paired SNR ratio 2.01 (PCCT is 2x better). PCCT noise 32% lower (mean 21.5 HU vs 32.8 HU). PCCT lower noise in 8/9 paired cases. N=28 PCCT, N=9 EID, N=9 paired.

---

### 1.4 Reconstruction Kernel

- **Threshold:** Soft-tissue kernel available; sharp kernel documented
- **Method:** Protocol review + DICOM header
- **Owner:** Murali
- **Status:** Complete
- **Evidence:** https://elucidbio.sharepoint.com/:w:/r/Technology/Documents/IP/PCCT/PT142_combined_analysis_COMPLETE_FINAL_v1.0.docx?d=w0fe5c7e3207c4790930cac6308e7bde7&csf=1&web=1&e=mQv35o
- **Notes:** Kernel family influences local vessel segment and patient level plaque characterization. Recon settings should be standardized for subsequent evaluations. Primary kernel: Bv44u\4 (soft-tissue, 120 kVp). All 30 PCCT cases are NAEOTOM Alpha, all EID cases are Somatom Force Bv40d\3.
