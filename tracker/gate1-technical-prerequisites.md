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

- **Threshold:** Peak aortic HU >= 250 in >= 90% of cases
- **Method:** Aortic ROI measurement (10mm cube at aortic centroid via compute_snr.py)
- **Owner:**
- **Status:** PCCT PASS / EID PASS
- **Evidence:** gate_results/snr_pcct.csv, gate_results/snr_eid.csv
- **Notes:** PCCT: 28/28 (100%) >= 250 HU -- PASS. EID: 24/25 (96%) >= 250 HU -- PASS. Only PT-124 (242 HU) below threshold. Threshold lowered from 300 to 250 HU on 2026-05-07 to reflect realistic scanner-side acceptance for diagnostic-quality CTA -- 250 HU corresponds to the lower bound of acceptable lumen-to-wall contrast for plaque characterization.

---

### 1.3 Image Noise (SNR/CNR)

- **Threshold:** Non-inferior to reference EID CTA (SNR ratio >= 0.85)
- **Method:** Aortic ROI SNR = mean_HU / std_HU (10mm cube at aortic centroid)
- **Owner:**
- **Status:** PASS
- **Evidence:** gate_results/snr_pcct.csv, gate_results/snr_eid.csv
- **Notes:** Paired SNR ratio 2.04 (PCCT ~2x better). PCCT lower noise in 24/25 paired cases. N=28 PCCT, N=25 EID. PT-142 EID Layne (wi-e5b7995e, wi-61c53cda), PT-152 EID, PT-158 EID still pending -- Aorta.nrrd not generated. PT-142 EID now uses Mackenzie's wi-ba36c00d (only ready PT-142 EID workitem).

---

### 1.4 Reconstruction Kernel

- **Threshold:** Soft-tissue kernel available; sharp kernel documented
- **Method:** Protocol review + DICOM header
- **Owner:** Murali
- **Status:** Complete
- **Evidence:** https://elucidbio.sharepoint.com/:w:/r/Technology/Documents/IP/PCCT/PT142_combined_analysis_COMPLETE_FINAL_v1.0.docx?d=w0fe5c7e3207c4790930cac6308e7bde7&csf=1&web=1&e=mQv35o
- **Notes:** Kernel family influences local vessel segment and patient level plaque characterization. Recon settings should be standardized for subsequent evaluations. Primary kernel: Bv44u\4 (soft-tissue, 120 kVp). All 30 PCCT cases are NAEOTOM Alpha, all EID cases are Somatom Force Bv40d\3.
