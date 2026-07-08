# Version 1 — Original analysis (FROZEN SNAPSHOT)

Frozen snapshot of the **original** PCCT qualification analysis, extracted from
git commit **`7fe43dc`** ("Add sub-segment intersection analysis as parallel Gate
3/4 pass"), the last state before the 2026-07 methodology corrections and data
refresh.

**Open `qualification_report_v1_original.html`** for the comprehensive original
report (self-contained, embedded figures + detailed results). *(Renamed from
`qualification_report.html` so it doesn't collide with the current v2 report of
the same name in `../gate_results/`.)*

## What this version is
| Property | Value |
|---|---|
| wCV estimator | **legacy `rms-rel`** — `sqrt(mean(d²/(2m²)))×100` (log branch = relative CV of logs; later found non-OQ-consistent) |
| Gate 4 bias criterion | project-specific `|bias| < 5%/10% of mean` |
| Input data | **original workitem summaries (pre-2026-07-07 re-work)** |
| N | 25 paired patients |
| Region | canonical named-vessel overlap + sub-segment parallel pass |

## ⚠️ Not re-runnable
`workitem_summaries/` was never version-controlled and the original CSVs were
overwritten by the 2026-07-07 refresh, so the original **input data is gone**.
This folder is a frozen record of the original *results, report, and figures*
only — it cannot be regenerated. The current (re-runnable) analysis lives in
`../gate_results/`.

## Contrast with Version 2 (current)
See `../tracker/statistical-methodology.md` for what changed (corrected
variance-component wCV, `--scanner-term`, OQ-consistent `--bias-criterion`,
verified component definitions) and `../gate_results/variants/` for the current
per-configuration results.
