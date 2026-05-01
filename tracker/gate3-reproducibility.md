# Gate 3 -- Core Quantitative Reproducibility

> **Performance** -- Scanner conditionally failed on miss; may expand N or adjust acquisition protocol and re-test once

Assessed at patient level using wCV (within-subject coefficient of variation) from paired PCCT/EID acquisitions of the same patient. Volumes are length-normalized before computing wCV to isolate variability due to image modality from vessel extent differences.

**Acceptance criterion:** 95% CI overlap between PCCT wCV and B.1P Delta Validation OQ wCV (730-CVV-040). Cross-scanner variability must be non-inferior to validated inter-operator variability.

**Primary metric:** Log-wCV for volume endpoints. Untransformed for length (descriptive only).

**Reference studies:**
- Original B.1P OQ (4-B1P-033 v2.0): log-transformed for process outputs, untransformed for plaque volumes
- Delta Validation OQ (730-CVV-040 v0.1): log-wCV with 95% bootstrap CIs

**Statistical design:** N >= 30 paired patients. wCV = sqrt(mean(d^2 / (2m^2))) x 100. 95% CI via 2000-sample bootstrap.

---

### 3.1 Lumen Volume

- **Reference:** Delta OQ log-wCV 7.5% [5.9%, 9.3%]; Original OQ 7.64% [5.52%, 8.55%]
- **Owner:**
- **Status:** PASS (N=9, preliminary)
- **Evidence:** gate_results/gate_summary.txt
- **Notes:** Log-wCV 12.97% [8.98%, 16.21%]. CI overlaps with delta OQ. Elevated point estimate expected for cross-scanner comparison vs inter-operator.

---

### 3.2 Wall Volume

- **Reference:** Delta OQ log-wCV 13.1% [10.2%, 16.0%]; Original OQ 9.51% [6.60%, 11.16%]
- **Owner:**
- **Status:** PASS (N=9, preliminary)
- **Evidence:** gate_results/gate_summary.txt
- **Notes:** Log-wCV 24.91% [13.78%, 35.14%]. CI overlaps with delta OQ. Wide CI at N=9 will narrow with more data. Wall collapse design change (730-CVV-017) applies equally to both scan types.

---

### 3.3 Vessel Volume

- **Reference:** Delta OQ log-wCV 8.7% [6.8%, 10.6%]; Original OQ 6.03% [4.44%, 6.74%]
- **Owner:**
- **Status:** PASS on log / FAIL on untransformed (N=9, preliminary)
- **Evidence:** gate_results/gate_summary.txt
- **Notes:** Log-wCV 12.44% [8.48%, 15.73%] -- CI overlaps delta OQ on log scale. Untransformed CI [20.29%, 36.66%] does not overlap [7.83%, 12.38%] -- consistent with vessel volume inheriting wall variability (VesselVol = Lumen + Wall).

---

### 3.4 Vessel Length (Descriptive)

- **Reference:** Delta OQ untransformed 13.88% [10.72%, 17.38%]; Original OQ 5.27% [3.43%, 6.52%]
- **Owner:**
- **Status:** Descriptive only (not acceptance criterion)
- **Evidence:** gate_results/gate_summary.txt
- **Notes:** Untransformed wCV 44.37% [28.94%, 55.46%]. Elevated due to different vessel tracing extent between PCCT and EID scans. Per 730-CVV-040 recommendation, length variability reflects analyst termination judgment, not segmentation error. Controlled by length normalization in volume analyses.

---

### 3.5 CALC Volume

- **Reference:** Delta OQ log-wCV 4.5% [3.5%, 5.4%]; Original OQ 13.9% (untransformed)
- **Owner:**
- **Status:** FAIL (N=9, preliminary)
- **Evidence:** gate_results/gate_summary.txt
- **Notes:** Log-wCV 26.17% [13.27%, 35.68%]. CI does not overlap delta OQ [3.5%, 5.4%]. CALC segmentation is HU-threshold sensitive and may be affected by kVp differences between PCCT (120 kVp) and EID (90-130 kVp).

---

### 3.6 LRNC Volume

- **Reference:** Delta OQ log-wCV 5.4% [4.1%, 6.9%]; Original OQ 59.3% (untransformed)
- **Owner:**
- **Status:** FAIL (N=9, preliminary)
- **Evidence:** gate_results/gate_summary.txt
- **Notes:** Log-wCV 80.52% [54.67%, 104.07%]. Highly elevated. LRNC delineation is the most challenging plaque component on CT and is sensitive to differences in image texture between PCCT and EID photon-counting vs energy-integrating detectors.

---

### 3.7 NonCALC Matrix Volume

- **Reference:** Delta OQ log-wCV 13.6% [10.8%, 16.4%]; Original OQ 58.6% (untransformed, Total NCP)
- **Owner:**
- **Status:** FAIL (N=9, preliminary)
- **Evidence:** gate_results/gate_summary.txt
- **Notes:** Log-wCV 35.22% [19.69%, 48.21%]. CI does not overlap delta OQ on log scale. May improve with more data.

---

### 3.8 Total Plaque Volume

- **Reference:** Delta OQ log-wCV 13.2% [10.4%, 16.0%]; Original OQ 44.1% (untransformed)
- **Owner:**
- **Status:** FAIL (N=9, preliminary)
- **Evidence:** gate_results/gate_summary.txt
- **Notes:** Log-wCV 27.37% [16.16%, 38.63%]. CI lower bound (16.16%) overlaps with delta OQ upper bound (16.0%) -- marginal. May pass with more paired data.
