# Case Review Notes

Per-patient observations from visual review, captured to inform interpretation of Gate 3/4 results. Distinguishes scanner-modality effects (PCCT vs EID) from acquisition-condition effects (motion, contrast timing, image quality).

---

### PT-133

- **Reviewer:**
- **Review date:** 2026-05-08
- **Findings:**
  - PCCT has motion artifact in the right coronary territory, leading to **missed plaque on the right side**. EID does not have motion artifact in that area.
  - Bad motion artifact on EID causes **overcall of wall** and **false NonCALC Matrix**. PCCT (no artifact) shows the plaque boundary more clearly.
- **Interpretation:** Discrepancy in Wall and NonCALC Matrix is a product of **different imaging time points and raw data**, not a PCCT-vs-EID modality difference. The motion artifact happens to fall on opposite vessels in the two acquisitions.
- **Implication for analysis:** Excluding PT-133 from sensitivity analysis is justified because the per-vessel artifact distribution is acquisition-specific, not scanner-specific.

---

### PT-142

- **Reviewer:**
- **Review date:** 2026-05-08
- **Findings:**
  - Wall is generally more expanded in EID.
  - Outer plaque boundary is much more defined in PCCT.
- **Interpretation:** Probably a **true bias at high disease burden** — when wall boundary is ambiguous, even slightly degraded image quality (EID's lower SNR / softer edges) shifts the segmentation outward. PCCT's higher SNR resolves the outer wall edge more reliably.
- **Implication for analysis:** Wall over-call on EID at high disease is likely a real modality-driven effect, not an acquisition artifact. Consistent with the persistent Wall log-wCV FAIL in Gate 3 and NonCALC Matrix proportional-bias FAIL in Gate 4.
