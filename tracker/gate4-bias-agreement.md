# Gate 4 -- Systematic Bias & Agreement

> **Performance** -- Flag for clinical review on failure; may restrict use to specific indications

**Method:** Bland-Altman on log(x+1) scale using raw (not length-normalized) volumes for direct comparability with delta OQ reference (730-CVV-040 Table 6).

**Proportional bias:** r^2 with 95% bootstrap CI. PASS if CI lower bound < 0.1.

---

### 4.1 Mean Bias -- Lumen Volume

- **Threshold:** |bias| < 5% of mean lumen volume
- **Method:** Bland-Altman plot on log(x+1) scale; report LoA
- **Owner:**
- **Status:** FAIL on bias / PASS on proportional bias (N=9, preliminary)
- **Evidence:** gate_results/bland_altman_plots/BA_LumenVol.png, gate_results/gate_summary.txt
- **Notes:** Log-scale bias 0.22 (PCCT > EID). Untransformed bias 333 mm3 (20.9% of mean) -- exceeds 5% threshold. Proportional bias r2=0.013 [0.000, 0.661] -- PASS. Bias driven by PCCT having longer vessel tracing (more total volume). This is a vessel extent difference, not a segmentation bias.

---

### 4.2 Mean Bias -- CALC Volume

- **Threshold:** |bias| < 10% of mean calc volume
- **Method:** Bland-Altman plot on log(x+1) scale; proportional bias test
- **Owner:**
- **Status:** FAIL on bias / PASS on proportional bias (N=9, preliminary)
- **Evidence:** gate_results/bland_altman_plots/BA_CALCVol.png
- **Notes:** Log-scale bias 0.52. Untransformed bias 58 mm3 (42.7% of mean). Proportional bias r2=0.090 [0.003, 0.567] -- PASS. Delta OQ ref bias: 0.02, ref LoA: [-0.17, 0.21]. PCCT systematically reports higher CALC -- potentially related to HU threshold sensitivity at different kVp.

---

### 4.3 Mean Bias -- LRNC Volume

- **Threshold:** |bias| < 10% of mean
- **Method:** Bland-Altman on log(x+1) scale
- **Owner:**
- **Status:** FAIL on bias / PASS on proportional bias (N=9, preliminary)
- **Evidence:** gate_results/bland_altman_plots/BA_LRNCVol.png
- **Notes:** Log-scale bias 0.86 with wide LoA [-1.45, 3.17]. Proportional bias r2=0.055 [0.000, 0.855] -- PASS. LRNC is inherently the most variable plaque component.

---

### 4.4 Mean Bias -- NonCALC Matrix Volume

- **Threshold:** |bias| < 10% of mean
- **Method:** Bland-Altman on log(x+1) scale
- **Owner:**
- **Status:** PASS on bias / PASS on proportional bias (N=9)
- **Evidence:** gate_results/bland_altman_plots/BA_NonCALCMATXVol.png
- **Notes:** Log-scale bias 0.08, LoA [-0.80, 0.97]. Untransformed bias 2.67 (1.0% of mean) -- PASS. Proportional bias r2=0.225 [0.019, 0.808] -- PASS. Delta OQ ref LoA: [-0.82, 0.67].

---

### 4.5 Mean Bias -- Total Plaque Volume

- **Threshold:** |bias| < 10% of mean
- **Method:** Bland-Altman on log(x+1) scale
- **Owner:**
- **Status:** FAIL on bias / PASS on proportional bias (N=9, preliminary)
- **Evidence:** gate_results/bland_altman_plots/BA_TotalPlaqueVolume.png
- **Notes:** Log-scale bias 0.25, LoA [-0.59, 1.09]. Untransformed bias 72.8 (17.6% of mean). Proportional bias r2=0.188 [0.005, 0.765] -- PASS. Delta OQ ref LoA: [-0.89, 0.78].

---

### 4.6 Limits of Agreement Summary

- **Threshold:** LoA within +/-1.96 SD; no proportional bias (r2 CI lower bound < 0.1)
- **Owner:**
- **Status:** All 7 variables PASS proportional bias test (N=9)
- **Evidence:** gate_results/gate_summary.txt
- **Notes:** All r2 95% CI lower bounds < 0.1. No evidence of proportional bias in any variable at current N. LoA will tighten with more paired data.
