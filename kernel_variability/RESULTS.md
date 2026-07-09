# Convolution-kernel + auto-segmentation variability — PT-142 PCCT

**Question.** For an *automatically* processed case (no reader edits), the plaque
variability can only come from (a) the **convolution kernel** and (b) the
**auto-segmentation algorithm**. How large are those, and what is the *pure kernel*
effect once the segmentation's variable traced extent is removed?

**Data.** One PT-142 PCCT acquisition reconstructed with **20 convolution kernels**
(families **Qr** ×10 and **Bv** ×10; sharpness 36/40/44; QIR strength 2/3/4; `u`/`d`
mode), each auto-segmented. Patient-level totals = sum of Left + Right coronary
targets per reconstruction. Source: `PT-142_workitem_summary_combined.xlsx`
(root; git-ignored). Reproduce with `python analyze_kernel_variability.py`.

> Single subject — this characterizes reconstruction/segmentation *sensitivity* for
> one patient, not a population estimate. LRNC ≈ 0 mm³ in this patient, so its metrics
> are near-zero noise and are excluded from interpretation.

---

## 1. Raw across-kernel variability (kernel + auto-seg combined)

| Endpoint (raw total) | mean | CV% | min–max | range/mean |
|---|--:|--:|--:|--:|
| Total Plaque | 966 mm³ | **2.4%** | 933–1010 | 8% |
| CALC | 667 mm³ | **3.0%** | 636–711 | 11% |
| Wall | 1123 mm³ | 3.5% | 1062–1186 | 11% |
| NonCALC Matrix | 299 mm³ | 9.4% | 253–361 | 36% |
| Lumen | 4153 mm³ | 9.8% | 2913–4490 | 38% |
| **Len (traced extent)** | 1224 mm | **21.3%** | 717–1688 | 79% |

Whole-vessel **totals are stable** (2–3% CV for total plaque and calcium), but the
**traced vessel length swings 79%** — the largest single source of variability, and a
property of the auto-segmentation, not of plaque content.

## 2. Separating extent from kernel — variance decomposition

Using the exact identity `log(V) = log(V/Len) + log(Len)`:

| Metric | SD(log V) *total* | SD(log Len) *extent* | SD(log density) *per-mm* | corr(Len, density) |
|---|--:|--:|--:|--:|
| Total Plaque | 2.4% | 21.9% | 22.2% | **−0.99** |
| CALC | 3.0% | 21.9% | 19.7% | **−1.00** |
| NonCALC Matrix | 9.4% | 21.9% | 28.9% | −0.97 |
| Wall | 3.5% | 21.9% | 22.6% | −0.99 |
| Lumen | 11.1% | 21.9% | 18.4% | −0.86 |

**Extent and per-mm density are almost perfectly anti-correlated (r ≈ −1).** When the
segmentation traces more vessel it adds low-plaque distal length, so the per-mm density
drops proportionally and the raw total barely moves. **The ~22% per-mm variability is a
mechanical artifact of the variable traced extent, not the kernel reclassifying tissue.**
Length-normalizing this data therefore *inflates* apparent variability.

## 3. Pure kernel effect — extent-invariant composition ratios

Dimensionless ratios (traced extent cancels) isolate what the kernel does to the
tissue *classification*:

| Ratio | mean | CV% | Qr − Bv |
|---|--:|--:|--:|
| CALC fraction of plaque | 0.690 | **3.5%** | −5% |
| CALC / Wall | 0.594 | 4.9% | −9% |
| Plaque burden `TP/(L+W)` | 0.185 | 10.3%\* | −3% |
| Wall / Lumen | 0.273 | 11.4%\* | +1% |
| **NonCALC fraction of plaque** | 0.309 | **8.0%** | **+12%** |
| NonCALC Matrix / Wall | 0.266 | 7.6% | +8% |

- **Calcium is kernel-robust** — the calcified fraction of plaque is the most stable
  quantity (3.5% CV); the kernel does not meaningfully change how much calcium is found.
- **Soft / non-calcified plaque is the one kernel-sensitive quantity** — ~8% CV, and
  **systematically ~+12% higher on the Qr family than Bv** (Qr assigns more tissue to
  soft matrix).

\* The plaque-burden and wall/lumen ratios keep a small residual extent coupling: the two
short-tracing `d_3` kernels (Len ≈ 717–831 mm) measure mostly the diseased proximal
segment, inflating those ratios. A *fully* extent-matched separation requires the
centerline-resolved data (a common distance-from-ostium window), not this patient-level
summary.

---

## Bottom line

For auto-processed PCCT of PT-142:

1. **Total plaque burden and calcium are reproducible across reconstruction kernels to
   ~2–3% CV.**
2. The dominant raw-volume variability (~22%) is the **auto-segmentation's traced extent**,
   which the kernel drives (smoother Bv traces ~23% longer than Qr) and which cancels out
   of the total.
3. The **pure convolution-kernel effect** on plaque content is small and specific: a
   modest, systematic shift of the **calcified↔soft split — soft-matrix fraction ~8% CV
   and ~+12% higher on Qr vs Bv** — with calcium essentially unchanged.

## Files

| File | Contents |
|---|---|
| `analyze_kernel_variability.py` | Reproducible analysis (reads the combined workbook) |
| `outputs/per_kernel_patient_level.csv` | Per-kernel patient-level totals, per-mm, ratios, parsed factors |
| `outputs/variability_summary.csv` | Raw + per-mm variability per endpoint |
| `outputs/variance_decomposition.csv` | log-identity split: extent vs density + correlation |
| `outputs/composition_ratios.csv` | Extent-invariant ratios: CV% + Qr-vs-Bv family effect |

## Next steps (not yet done)

- **True fixed-extent kernel effect** via PT-142's per-kernel centerline/NRRD outputs
  (common distance-from-ostium window), removing the residual proximal/distal coupling.
- **Extend to more patients** with multi-kernel reconstructions for a population estimate.
