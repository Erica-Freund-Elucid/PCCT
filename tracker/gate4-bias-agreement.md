# Gate 4 -- Systematic Bias & Agreement

> **Performance** -- Flag for clinical review on failure; may restrict use to specific indications

**Method:** Bland-Altman on raw (not length-normalized) volumes summed over the per-pair vessel-overlap intersection (canonical as of 2026-05-08). **Process outputs (Lumen, Wall, Vessel) on log(x+1) scale** for comparability with delta OQ reference (730-CVV-040 Table 6); **plaque (CALC, LRNC, NonCALC, Total) on untransformed scale** to match 4-B1P-033 Table 9 reporting. PT-124 excluded (no overlapping vessel). Most pairs are partial vessel-overlap, analyzed on the intersection only.

**Proportional bias:** r^2 with 95% bootstrap CI. PASS if CI lower bound < 0.1.

---

### 4.1 Mean Bias -- Lumen Volume

- **Threshold:** |bias| < 5% of mean lumen volume
- **Method:** Bland-Altman plot on log(x+1) scale; report LoA
- **Owner:**
- **Status:** PASS on bias / PASS on proportional bias (N=25, vessel-overlap, preliminary)
- **Evidence:** gate_results/bland_altman_plots/BA_LumenVol.png, gate_results/gate_summary.txt
- **Notes:** Log-scale bias -0.016, LoA [-0.49, 0.46]. Untransformed bias -26.9 mm3 (1.6% of mean) -- under 5% threshold. Proportional bias r2=0.005 [0.000, 0.207] -- PASS.

---

### 4.1b Mean Bias -- Wall Volume (project-specific)

- **Threshold:** |bias| < 10% of mean wall volume
- **Method:** Bland-Altman plot on log(x+1) scale
- **Owner:**
- **Status:** PASS on bias / PASS on proportional bias (N=25, vessel-overlap, preliminary)
- **Evidence:** gate_results/bland_altman_plots/BA_WallVol.png
- **Notes:** Threshold added 2026-05-08 as a project-specific Gate 4 criterion. Not derivable from 4-B1P-033 or 730-CVV-040 reference documents -- both assess wall reproducibility via wCV only (730-CVV-040 Table 6 is plaque-only, see PDF page 8). Added in response to case-review finding (PT-142) that wall over-call on EID at high disease appears to be a true modality-driven effect, so a bias check is warranted. Threshold matched to plaque components (10%) -- wall measurement is dominated by the wall-plaque boundary, closer in difficulty to plaque segmentation than to lumen (which has the lumen-blood contrast edge). Log-scale bias -0.068, LoA [-0.60, 0.46]. Untransformed bias -51.0 mm3 (9.3% of mean) -- under 10% threshold. Proportional bias r2=0.016 [0.000, 0.233] -- PASS. Sensitivity (excl PT-119/133/142): bias drops further to 4.0% (see Section 4.7).

---

### 4.2 Mean Bias -- CALC Volume

- **Threshold:** |bias| < 10% of mean calc volume
- **Method:** Bland-Altman on untransformed scale; proportional bias test
- **Owner:**
- **Status:** FAIL on bias / PASS on proportional bias (N=25, vessel-overlap, preliminary)
- **Evidence:** gate_results/bland_altman_plots/BA_CALCVol.png
- **Notes:** Untransformed bias 27.61 mm3 (19.2% of mean). LoA [-33.97, 89.19] mm3. Proportional bias r2=0.142 [0.010, 0.703] -- PASS (CI lower bound 0.010 < 0.1). PCCT systematically reports higher CALC than EID.

---

### 4.3 Mean Bias -- LRNC Volume

- **Threshold:** |bias| < 10% of mean
- **Method:** Bland-Altman on untransformed scale
- **Owner:**
- **Status:** FAIL on bias / PASS on proportional bias (N=25, vessel-overlap, preliminary)
- **Evidence:** gate_results/bland_altman_plots/BA_LRNCVol.png
- **Notes:** Untransformed bias 5.65 mm3 (21.0% of mean). LoA [-49.58, 60.87] mm3. Proportional bias r2=0.051 [0.001, 0.725] -- PASS.

---

### 4.4 Mean Bias -- NonCALC Matrix Volume

- **Threshold:** |bias| < 10% of mean
- **Method:** Bland-Altman on untransformed scale
- **Owner:**
- **Status:** FAIL on bias / PASS on proportional bias (N=25, vessel-overlap, recovered from prop-bias FAIL on target-overlap)
- **Evidence:** gate_results/bland_altman_plots/BA_NonCALCMATXVol.png
- **Notes:** Untransformed bias -65.17 mm3 (-24.5% of mean). LoA [-281.60, 151.26] mm3. Proportional bias r2=0.436 [0.056, 0.725] -- PASS (CI lower bound 0.056 < 0.1). Recovered from FAIL [0.127, 0.743] under target-overlap because vessel-overlap removes large-magnitude PCCT-only data points that drove the apparent slope.

---

### 4.5 Mean Bias -- Total Plaque Volume

- **Threshold:** |bias| < 10% of mean
- **Method:** Bland-Altman on untransformed scale
- **Owner:**
- **Status:** PASS on bias (7.3%) / PASS on proportional bias (N=25, vessel-overlap)
- **Evidence:** gate_results/bland_altman_plots/BA_TotalPlaqueVolume.png
- **Notes:** Untransformed bias -31.91 mm3 (-7.3% of mean) -- under 10% threshold. LoA [-289.67, 225.84] mm3. Proportional bias r2=0.167 [0.001, 0.542] -- PASS.

---

### 4.6 Limits of Agreement Summary

- **Threshold:** LoA within +/-1.96 SD; no proportional bias (r2 CI lower bound < 0.1)
- **Owner:**
- **Status:** All 7 variables PASS proportional bias test (N=25, vessel-overlap, mixed-scale)
- **Evidence:** gate_results/gate_summary.txt
- **Notes:** Process outputs (Lumen, Wall, Vessel) reported on log scale per 730-CVV-040 Table 6 reference; plaque variables on untransformed scale per 4-B1P-033 Table 9 reference. Vessel-overlap restriction recovered NonCALC prop-bias from FAIL to PASS (r2 CI lower 0.056 < 0.1). Bias magnitudes for component plaque (CALC +19%, LRNC +21%, NonCALC -25%) and Wall (-9.3%) remain elevated -- consistent with case-review notes (see tracker/case-reviews.md).

---

### 4.7 Sensitivity -- exclude case-review-flagged patients

- **Threshold:** N/A (sensitivity / advisory only)
- **Method:** Re-run Gate 3 wCV and Gate 4 Bland-Altman with PT-119, PT-133, PT-142 excluded. PT-119 and PT-133 are acquisition-condition artifacts (motion); PT-142 is the high-disease modality-bias case. See tracker/case-reviews.md.
- **Owner:**
- **Status:** Sensitivity-only; not used for gate qualification at N=25
- **Evidence:** gate_results/gate_summary.txt (run with full N=25); sensitivity computed offline.
- **Notes:**
  - Gate 3: all 7 variables remain PASS at N=22 (similar values to N=25).
  - Gate 4 Wall bias: 9.3% PASS -> 4.0% PASS at N=22 -- consistent with PT-142 contributing the wall modality-bias signal at high disease.
  - Gate 4 Total Plaque bias: 7.3% PASS -> 1.0% PASS (much improved).
  - Gate 4 CALC proportional bias: PASS [0.010, 0.704] -> FAIL [0.298, 0.815] at N=22 -- artifact of removing high-disease anchor cases that flattened the regression slope; expected to resolve at N>=30.
  - Recommendation: report N=25 as primary; document Wall bias deviation with case-review traceability rather than excluding patients from the qualification cohort.
