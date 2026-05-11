# Gate 3 -- Core Quantitative Reproducibility

> **Performance** -- Scanner conditionally failed on miss; may expand N or adjust acquisition protocol and re-test once

Assessed at patient level using wCV (within-subject coefficient of variation) from paired PCCT/EID acquisitions of the same patient. Volumes are length-normalized before computing wCV to isolate variability due to image modality from vessel extent differences.

**Acceptance criterion:** 95% CI overlap between PCCT wCV and B.1P Delta Validation OQ wCV (730-CVV-040). Cross-scanner variability must be non-inferior to validated inter-operator variability.

**Metric per endpoint:** Log-wCV for process outputs (Lumen, Wall, Vessel vol -- matches 4-B1P-033 Table 7). Untransformed wCV for plaque volumes (matches 4-B1P-033 Table 9 reporting). Untransformed wCV for length (descriptive only).

**Reference studies:**
- Original B.1P OQ (4-B1P-033 v2.0): log-transformed for process outputs, untransformed for plaque volumes
- Delta Validation OQ (730-CVV-040 v0.1): log-wCV with 95% bootstrap CIs

**Statistical design:** N >= 30 paired patients. wCV = sqrt(mean(d^2 / (2m^2))) x 100. 95% CI via 2000-sample bootstrap.

**Per-pair vessel-overlap restriction (canonical as of 2026-05-08):** for each PCCT/EID pair, only the (bodySite, vessel-location) intersection is summed -- so the same anatomical extent is measured on each scanner. Pairs with no overlapping vessel are excluded entirely. PT-124 excluded (PCCT=Left vessels, EID=Right vessel only). Most pairs are partial-overlap; PCCT typically traces 1-3 additional distal vessels (PDA, PLB, Diagonals, marginals). Restricting to the intersection eliminates Length-normalization dilution from PCCT's longer trace -- this is what flipped Wall log-wCV from FAIL (target-overlap) to PASS at N=25.

---

### 3.1 Lumen Volume

- **Reference:** Delta OQ log-wCV 7.5% [5.9%, 9.3%]; Original OQ 7.64% [5.52%, 8.55%]
- **Owner:**
- **Status:** PASS (N=25, vessel-overlap, preliminary)
- **Evidence:** gate_results/gate_summary.txt
- **Notes:** Log-wCV 8.72% [6.99%, 10.30%]. CI overlaps delta OQ [5.9%, 9.3%]. Tightened from target-overlap [8.27%, 12.28%].

---

### 3.2 Wall Volume

- **Reference:** Delta OQ log-wCV 13.1% [10.2%, 16.0%]; Original OQ 9.51% [6.60%, 11.16%]
- **Owner:**
- **Status:** PASS (N=25, vessel-overlap, preliminary -- recovered from FAIL on target-overlap)
- **Evidence:** gate_results/gate_summary.txt
- **Notes:** Log-wCV 18.48% [14.07%, 22.55%]. CI lower bound (14.07%) below delta OQ upper bound (16.0%) -- overlaps. Was 23.58% [17.18%, 30.02%] FAIL on target-overlap; switching to vessel-overlap eliminated the artifactual length-normalization dilution caused by PCCT tracing additional thin-walled distal vessels.

---

### 3.3 Vessel Volume

- **Reference:** Delta OQ log-wCV 8.7% [6.8%, 10.6%]; Original OQ 6.03% [4.44%, 6.74%]
- **Owner:**
- **Status:** PASS on log / FAIL on untransformed (N=25, vessel-overlap, preliminary)
- **Evidence:** gate_results/gate_summary.txt
- **Notes:** Log-wCV 8.40% [6.83%, 9.85%] -- CI overlaps delta OQ on log scale. Untransformed CI [16.67%, 24.24%] does not overlap [7.83%, 12.38%] -- consistent with raw-scale heteroscedasticity at high disease.

---

### 3.4 Vessel Length (Descriptive)

- **Reference:** Delta OQ untransformed 13.88% [10.72%, 17.38%]; Original OQ 5.27% [3.43%, 6.52%]
- **Owner:**
- **Status:** Descriptive only (not acceptance criterion)
- **Evidence:** gate_results/gate_summary.txt
- **Notes:** Untransformed wCV 24.64% [18.88%, 29.60%] (down from 33.6% with target-overlap). Mean PCCT-EID length difference dropped from +44% to +26% on overlapping vessels -- residual difference reflects analyst termination judgment within shared vessels.

---

### 3.5 CALC Volume

- **Reference:** Delta OQ untransformed wCV 25.78% [18.70%, 35.75%]; Original OQ 13.9% (untransformed)
- **Owner:**
- **Status:** PASS (N=25, vessel-overlap, preliminary)
- **Evidence:** gate_results/gate_summary.txt
- **Notes:** Untransformed wCV 21.08% [16.57%, 25.03%]. CI overlaps delta OQ [18.70%, 35.75%]. Tightened from target-overlap [18.65%, 30.91%].

---

### 3.6 LRNC Volume

- **Reference:** Delta OQ untransformed wCV 78.80% [48.09%, 163.40%]; Original OQ 59.3% (untransformed)
- **Owner:**
- **Status:** PASS (N=25, vessel-overlap, preliminary)
- **Evidence:** gate_results/gate_summary.txt
- **Notes:** Untransformed wCV 94.24% [78.29%, 108.83%]. CI overlaps delta OQ [48.09%, 163.40%]. LRNC inherently the most variable plaque component.

---

### 3.7 NonCALC Matrix Volume

- **Reference:** Delta OQ untransformed wCV 32.61% [24.87%, 41.94%]; Original OQ 58.6% (Total NCP, untransformed)
- **Owner:**
- **Status:** PASS (N=25, vessel-overlap, preliminary)
- **Evidence:** gate_results/gate_summary.txt
- **Notes:** Untransformed wCV 39.46% [30.22%, 47.73%]. CI lower bound (30.22%) below delta OQ upper bound (41.94%) -- overlaps.

---

### 3.8 Total Plaque Volume

- **Reference:** Delta OQ untransformed wCV 27.08% [20.68%, 34.74%]; Original OQ 44.1% (untransformed)
- **Owner:**
- **Status:** PASS (N=25, vessel-overlap, preliminary)
- **Evidence:** gate_results/gate_summary.txt
- **Notes:** Untransformed wCV 31.40% [24.65%, 37.47%]. CI lower bound (24.65%) below delta OQ upper bound (34.74%) -- overlaps. Tightened from target-overlap [28.58%, 46.69%].
