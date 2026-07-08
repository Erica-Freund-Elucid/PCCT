# Statistical Methodology & Verified Definitions

Consolidated methodology notes and empirically verified definitions for the PCCT
qualification analysis. Established/updated 2026-07-08.

Companion to `../run_gate_analyses.py` (implementation) and the per-gate trackers.

---

## 1. wCV estimator (Gate 3)

### What the OQ actually does (730-CVV-040 §11)
> "Statistical analysis was performed using R 4.3.2. Inter-run reproducibility
> was assessed using the within-subject coefficient of variation (%wCV)…
> **Log(x+1) transformation was applied to all variables prior to linear mixed
> model analysis**… Volumes were normalized to target length."

Cited method: **Quan & Shih (1996), "Assessing reproducibility by the
within-subject coefficient of variation with random effects models"** (Biometrics
52:1195). So the reference wCV is a **variance-component / random-effects (linear
mixed model)** estimator on **log(x+1)**, length-normalized volumes — not a
per-pair relative-CV formula.

### Correct estimator (implemented as `--wcv-method variance-component`, default)
For paired (k=2) data the within-subject variance component is `σ²_w = mean(d²/2)`
(d = PCCT − EID after length-normalization):
- **untransformed:** `wCV = sqrt(σ²_w) / grand_mean × 100`
- **log(x+1):** `wCV = sqrt(exp(σ²_w,log) − 1) × 100`

### Legacy estimator (`--wcv-method rms-rel`, retained for reproducing pre-2026-07 results)
`wCV = sqrt(mean(d²/(2·m²))) × 100`, m = (PCCT+EID)/2, applied to raw or log(x+1).
- **Untransformed** legacy ≈ variance-component (agrees within a few points).
- **Log branch is wrong:** it computes the *relative CV of the log-transformed
  values* (divides the log-difference by the **mean of the logs**), which is not a
  variance-component estimator and has no standard interpretation. It
  **understates** large-magnitude endpoints and **overstates** small ones.

### Impact of the fix (N=25, canonical named-vessel-overlap, log-wCV)
| Endpoint | legacy log-wCV | correct log-wCV | delta-OQ ref | verdict change |
|---|--:|--:|--:|---|
| Lumen | 7.2% | 14.1% | 7.5% [5.9, 9.3] | PASS → **FAIL** |
| Vessel | 8.4% | 18.5% | 8.7% [6.8, 10.6] | PASS → **FAIL** |
| Wall | 25.9% | 27.5% | 13.1% [10.2, 16.0] | FAIL → FAIL |

The comfortable Lumen/Vessel passes under the legacy formula were an artifact of
the non-standard log calculation.

---

## 2. Scanner term (`--scanner-term`, DEFAULT ON)

The PCCT−EID difference decomposes as
**scanner bias + scanner random dispersion + operator differences + traced-extent
differences + noise.** The OQ's wCV is *same-scanner inter-operator* dispersion
with **no systematic bias**. To compare like-for-like, remove the systematic
modality bias and compare only the random dispersion.

**This is the default and the standard basis for the Gate 3 wCV-vs-OQ acceptance:**
the OQ variability limit is a bias-free inter-operator dispersion, so the PCCT
statistic compared against it must also exclude the systematic scanner bias (which
is assessed separately in Gate 4). Use `--no-scanner-term` to see the raw,
bias-inflated cross-scanner wCV.

- **Implementation:** `σ²_w = Var(d)/2` (variance about the mean difference) instead
  of `mean(d²)/2`. The removed systematic bias is assessed separately in Gate 4.
- **Impact (N=25, canonical, log-wCV):** Lumen 14.2 → **11.2**, Wall 27.5 → **17.9**,
  Vessel 18.5 → **12.6** — all now marginally CI-overlap **PASS**. The process-output
  "failures" were driven mostly by systematic bias, not excess random dispersion.
- **Caveats:** with one read per patient·scanner (and different analysts / traced
  extents per scanner) the residual still confounds scanner with operator and
  extent → it is an **upper bound** on scanner-only variability, not a clean
  isolation. True isolation needs replicate reads (same patient+scanner, multiple
  operators). The sub-segment pass removes the extent confound.

---

## 2.1 Gate 4 bias acceptance criterion (`--bias-criterion`)

The systematic bias removed by the scanner term (§2) is exactly what Gate 4
should assess. Two criteria are available:

- **`pct-threshold`** (default, legacy): `|untransformed bias| < 5%` (lumen) /
  `10%` (wall, plaque) of the mean. **Project-specific — NOT derived from the OQ**
  (the Gate 4 tracker §4.1b admits this for wall; it holds for all).
- **`oq-ci-overlap`** (OQ-consistent): for the plaque endpoints, compute the PCCT
  **log(x+1)** Bland-Altman bias with a bootstrap **95% CI** and require it to
  **overlap the 730-CVV-040 Table 6 bias 95% CI** — the same CI-overlap philosophy
  as the Gate 3 wCV acceptance, on the OQ's own (log) scale. The OQ Table 6
  inter-operator bias CIs essentially include 0, so this asks "is the cross-scanner
  bias distinguishable from OQ inter-operator bias?". Process outputs have no OQ
  BA-bias reference (bias is *supporting*; wCV is primary per the OQ), so they are
  reported descriptively under this criterion.

OQ Table 6 references (log scale) wired into `GATE4_VARIABLES`:
CALC bias [−0.01, 0.05]; LRNC [−0.06, −0.01]; NonCALCMATX [−0.19, 0.05];
TotalPlaque [−0.19, 0.08].

**Result on latest data (N=25, canonical):**
| Endpoint | pct-threshold | oq-ci-overlap | note |
|---|---|---|---|
| CALC | PASS (2.7%) | PASS | bias ≈ 0 |
| LRNC | FAIL (17.2%) | PASS | CI huge [−0.70, 0.39] — **low power**, can't reject |
| NonCALC Matrix | FAIL (47%) | **FAIL** | PCCT lower, CI entirely below OQ |
| Total Plaque | FAIL (31%) | **FAIL** | PCCT lower, CI entirely below OQ |

**NonCALC Matrix and Total Plaque bias fail under both criteria — a real,
statistically distinguishable modality bias (PCCT systematically lower on the log
scale). This is the binding constraint.** LRNC "passes" the CI-overlap test only
because it is too noisy to reject (low power) — do not read that as agreement.
Wall bias is only tested under `pct-threshold` (project-specific, no OQ ref) and
has been growing across re-processing: 9.3% (May) → 20.7% (07-06) → 28.3% (07-07).

**Plaque Bland-Altman is on the log(x+1) scale** (changed 2026-07-08): the delta
OQ Table 6 plaque BA is log-scale, and raw plaque volumes are heteroscedastic
(spread grows with magnitude → funnel), so log is both OQ-consistent and the
statistically correct representation. The BA plots overlay the OQ Table 6 **bias
line, bias 95% CI band, and LoA band** against the **PCCT bias 95% CI band**, so
the CI-overlap acceptance is visible (title annotates PASS/FAIL). OQ-vs-PCCT
side-by-side tables with CIs are emitted to `gate_results/gate3_comparison.csv`
(wCV) and `gate4_comparison.csv` (bias) and rendered in the report.

Caveat: sub-segment bias verdicts are on **stale** segmentations; regenerate
before use.

---

## 3. Analysis region: canonical vs sub-segment

- **Canonical:** pairs by named-vessel `(bodySite, location)` overlap; whole-vessel
  totals, length-normalized. Within a shared vessel the two scanners may trace
  different arcs → residual extent confound.
- **Sub-segment intersection:** matches the two centerlines by distance-from-ostium
  (sub-voxel SDF) and recomputes volumes over the common arc only → identical
  anatomy on both scanners, extent confound removed. Runs as a parallel Gate 3/4
  pass on `workitem_summaries/subsegment/`. See `../scripts/subsegment/`.

**Stacking of the two controls (process outputs, log-wCV, N=25):**
| Variant | Lumen | Wall | Vessel | process outputs |
|---|--:|--:|--:|---|
| canonical, no scanner term | 14.2 | 27.5 | 18.5 | FAIL |
| canonical, + scanner term | 11.2 | 17.9 | 12.6 | PASS (marginal) |
| sub-segment, no scanner term | 11.7 | 12.5 | 10.6 | PASS |
| sub-segment, + scanner term | 10.0 | 12.0 | 9.1 | PASS (tightest) |

Sub-segment removes extent; scanner-term removes bias; together they give the
cleanest, all-passing process-output reproducibility.

---

## 4. Verified volume-component definitions (2026-07-08, from NRRD region labels)

Confirmed by reading the multi-NRRD ASCII headers on the ip3-manager instance.

**Plaque — `composition.multi.nrrd` has 7 region slices:**
`CALC, LRNC, IPH, PVAT, MATX, FIBL, NonCALCMATX`.
- `NonCALCMATX` is a **precomputed composite region (slice 6)**, NOT a sum of
  primitives (LRNC+IPH+MATX is ~34% off). Extract it directly.
- CSV columns: `CALCVol / LRNCVol / IPHVol / MATXVol / PVATVol` = those regions;
  `NonCALCMATXVol` = the NonCALCMATX region; `TotalNonCALCVolume` = NonCALCMATX.
- **`TotalPlaqueVolume = CALC + NonCALCMATX`** (≈3% SDF-merge/voxelization residual,
  not a different formula).

**Geometry — partition regions:**
- `lumenPartition.multi.nrrd` region = **`Lumen`**
- `wallPartition.multi.nrrd` region = **`LumenAndWall`** (the outer vessel boundary =
  lumen+wall, NOT wall-only) — confirmed across all analyst/date versions.

**CSV geometry columns** (identity `LumenAndWallVol = LumenVol + WallVol` is exact,
55/55 patients):
- `LumenVol` = lumen
- `LumenAndWallVol` = lumen + wall = **vessel**
- `WallVol` = `LumenAndWallVol − LumenVol` = wall-only tissue

**Sub-segment CSV builder mapping (reproduces the CSV):**
```
LumenVol        = pipeline raw "lumen"          (lumenPartition region = Lumen)
LumenAndWallVol = pipeline raw "wall"           (wallPartition region  = LumenAndWall)
WallVol         = raw "wall" − raw "lumen"       (wall-only tissue)
CALCVol         = CALC region        (direct)
NonCALCMATXVol  = NonCALCMATX region (direct)
TotalPlaqueVolume = CALC + NonCALCMATX
```
Because the pipeline's raw "wall" is the `LumenAndWall` region, `WallVol` must be
derived by subtracting lumen — mapping `WallVol = raw wall` would over-count by one
lumen and would NOT match the CSV. This confirms the builder mapping noted in
`../scripts/subsegment/RESUME.md`.

---

## 5. Open items / caveats

- **N=25, preliminary** (target ≥ 30). All process-output passes above are
  CI-overlap-driven with wide bootstrap CIs; point estimates still exceed the OQ
  point estimates.
- **Plaque bias (Gate 4) is real and within-vessel.** The sub-segment mitigation
  made plaque bias *worse*, showing the bias is concentrated within the shared
  centerline, not at the distal tail (see `../scripts/subsegment/RESUME.md`).
  The scanner term relocates bias to Gate 4; it does not remove it. Wall and Total
  Plaque bias still FAIL Gate 4.
- **Sub-segment CSVs are stale** relative to the 2026-07-07 workitem refresh
  (they reflect earlier segmentations). Regenerate via the sub-segment pipeline to
  reflect current data before citing sub-segment numbers as final.
- **Reference-threshold basis.** For a fully valid CI-overlap test the delta-OQ
  reference wCVs should be recomputed on the same variance-component basis; the OQ
  per-case data (Attachments 2 & 3 of 730-CVV-040) would be needed to confirm.
