# PCCT Sub-segment Intersection Mitigation — Resume Notes

**Last update:** 2026-05-22 (post-rescue). N=25 sub-segment data now matches canonical cohort (PT-124 excluded both, per tracker). **Mitigation did NOT work** even with full cohort — see "Findings" section.

## Goal

Mitigate failing Gate 4 Bland–Altman bias tests in `gate_results/gate_summary.txt`:

| Variable | Current bias | Threshold | State |
|---|---|---|---|
| CALC | **+19.2%** | \|bias\|<10% | FAIL |
| LRNC | **+21.0%** | \|bias\|<10% | FAIL |
| NonCALC Matrix | **−24.6%** | \|bias\|<10% | FAIL |

by recomputing volumes over **sub-segment intersection** of PCCT and EID centerlines
(distance-from-ostium match with subvoxel precision) — not just the current
(bodySite, location) vessel overlap. Source repo:
https://github.com/ElucidBioimaging/PCCTvsEID/tree/main/ComparePCCTandEIDWorkItems

**Per Erica:** Recompute **both** Gate 3 (wCV) and Gate 4 (BA) in parallel; output
goes to `workitem_summaries/subsegment/`; canonical sections stay untouched.

## What is done

- [x] Cloned mitigation scripts to `PCCT/scripts/subsegment/`:
  - `common_vessel_sections_plaque_volume.py` (62 KB)
  - `multi_nrrd_reader_writer.py` (17 KB)
  - `workitem_plaque_intersection_calculator.sh` (6 KB) — wrapper
  - `README.md` from upstream repo
- [x] Located workitem dirs on the dev box for PT-129:
  - PCCT (`wi-aff6fdcf`): `/inst/zenith/AppData/Working Storage/wi_layne.cassidy/PT-129/wi-aff6fdcf` (also under `wilist_`)
  - EID  (`wi-8622f0aa`): `/inst/zenith/AppData/Working Storage/wi_layne.cassidy/PT-129/wi-8622f0aa` (also under `wilist_`)
- [x] Verified dev box constraints:
  - `/tmp` is `tmpfs` **noexec** — venv must NOT live there (SimpleITK .so fails to load)
  - `/var/tmp` is on root fs (123 GB / 110 GB free), exec-able — use this
  - System python3 is PEP-668 externally-managed; must use venv
  - `venv` module is available
- [x] Created `/var/tmp/pcct_subseg/venv` with `numpy scipy SimpleITK pynrrd vtk` installed (last successful run before pause).
- [x] Cloned `PCCTvsEID` repo to `/var/tmp/pcct_subseg/PCCTvsEID/` (last successful run).

## What is NOT done

- [x] Verify the pipeline on PT-129 against pptx anchors — **PASSED 2026-05-22**.
  All 8 anchors matched exactly (Left CALC 180.64/131.81, Right CALC 82.30/43.18,
  Total 262.94/174.99, Diff 87.95, Mean 218.965). Wrapper wall time ~31s on the box.
  VTK `vtkMassProperties: Input data type must be VTK_TRIANGLE not 3` warnings are
  noisy but benign — the numbers come out right.
- [x] Replaced the failed `git clone https://github.com/ElucidBioimaging/PCCTvsEID`
  (repo requires auth on the box) with a tar+base64 upload of the local
  `PCCT/scripts/subsegment/` scripts. Payload file:
  `C:/Users/EricaFreund/AppData/Local/Temp/ssm_upload_scripts.json` (~27 KB).
  The wrapper script on the box is `workitem_plaque_intersection_calculator.sh`
  (not `run_calc_volume_pipeline.sh` as documented in the upstream README).
- [x] Run pipeline on the other 24 patients — 20 of 24 produced output (17 full,
  3 right-only). 4 left-only partials failed because the wrapper runs right first
  with `set -e`, so right-target absence kills the left pass too. PT-125 left also
  failed (duplicate `LeftCoronary` entries in its workitem.json).
- [x] Synced JSON outputs to `PCCT/workitem_summaries/subsegment/{PT-*}/final/{left,right}/...`.
- [x] Built per-patient sub-segment CSVs in
  `PCCT/workitem_summaries/subsegment/{PCCT,EID}/PT-*_workitem_summary_*_subseg.csv`.
  **WARNING**: sub-segment JSON's `wall` component is actually `LumenAndWall`
  (the wallPartition multi-NRRD slice region is `LumenAndWall`, not wall-only).
  CSV builder maps: `LumenAndWallVol = wall`, `WallVol = wall - lumen`,
  `LumenVol = lumen`. Confirmed by inspecting wallPartition NRRD slice metadata.
  **[VERIFIED 2026-07-08]** Re-confirmed against live NRRD headers on ip3-manager:
  `wallPartition` region = `LumenAndWall`, `lumenPartition` region = `Lumen`, and
  the CSV identity `LumenAndWallVol = LumenVol + WallVol` is exact (55/55 patients).
  This builder mapping is CORRECT — do NOT change to `WallVol = raw wall` (that
  would over-count by one lumen). Plaque: `composition.multi.nrrd` has regions
  `CALC, LRNC, IPH, PVAT, MATX, FIBL, NonCALCMATX`; `NonCALCMATX` is a precomputed
  composite slice (extract directly, do NOT sum primitives), and
  `TotalPlaqueVolume = CALC + NonCALCMATX`. Full write-up:
  `PCCT/tracker/statistical-methodology.md` §4.
- [x] Extended `run_gate_analyses.py` with a parallel sub-segment pass:
  `load_paired_data()` now accepts `pcct_dir, eid_dir` params; `run_gate4()` accepts
  `plot_dir`. Sub-segment Gate 3 + Gate 4 sections append to `gate_results/gate_summary.txt`.
  Sub-segment Bland-Altman plots go under `gate_results/bland_altman_plots_subsegment/`.

## Findings (2026-05-22, post-rescue N=25)

**Mitigation strategy did NOT work.** Matched-cohort comparison
(`gate_results/subsegment_comparison.txt`) on the same 25 patients
(canonical Gate 4 cohort):

| Variable | Canon bias / % | Sub-segment bias / % | Δ% | Result |
|---|---|---|---|---|
| CALC | 27.66 mm³ / 19.2% | 24.93 mm³ / 21.3% | +2.1 | **WORSE** |
| LRNC | 5.65 mm³ / 21.0% | 6.46 mm³ / 29.7% | +8.7 | **WORSE** |
| NCMATX | -65.20 mm³ / 24.6% | -76.62 mm³ / 36.6% | +12.1 | **WORSE** |
| TotalPlaque | -31.90 mm³ / 7.3% | -45.22 mm³ / 13.0% | +5.7 | **WORSE** |
| Lumen | -26.89 mm³ / 1.6% | -87.34 mm³ / 7.6% | +6.0 | WORSE (still <10%) |
| Wall (wall-only) | -51.05 / 9.3% | -42.07 / 9.2% | -0.1 | unchanged |
| LumenAndWall | -77.94 / 3.6% | -129.41 / 8.1% | +4.5 | WORSE |

PT-129 (the case-review pptx anchor) was the only patient with a large
reduction in CALC diff (119.89 → 87.95 mm³). PT-119 went the other way:
+52 → +104 mm³. PT-125 and PT-159 SIGN-FLIPPED their CALC bias under
sub-segment (PCCT and EID swap which is larger) — suggesting the
named-vessel-overlap measurement was being fooled by tail differences
that masked a real underlying bias.

The CALC bias absolute mm³ value barely changed (27.66 → 24.93).
The bias % went **up** because the means of both PCCT and EID dropped
more than the gap between them did, so |bias|/mean got worse.

**This invalidates the original hypothesis** that the bias is driven
by PCCT analysts tracing slightly farther into the same named vessel.
The bias is concentrated WITHIN the shared centerline, not at the
distal tail.

## Next steps (open)

## Stdout cap workaround for pulling JSON outputs

`aws ssm get-command-invocation` returns **at most 24000 chars** in
`StandardOutputContent`. The left `common_vessel_sections_summary.json` for
PT-129 is ~51 KB raw, so a plain `cat` truncates mid-file. Use gzip+base64
inline and decode locally:

```bash
gzip -c <file> | base64 -w0       # on the box; total b64 ≈ 5 KB for left
# locally:
python -c "import gzip,base64,json,sys; print(json.dumps(json.loads(gzip.decompress(base64.b64decode(sys.stdin.read().strip()))), indent=2))"
```

## Verification anchors for PT-129

From `PCCT/PCCT and CCTA intersection PT-129.pptx` (slide 7) — the script run
against PT-129 must match these to be accepted:

| Target | PCCT (wi-aff6fdcf) | EID (wi-8622f0aa) |
|---|---|---|
| Left CALC volume (mm³) | **180.64** | **131.81** |
| Right CALC volume (mm³) | **82.30** | **43.18** |
| Total CALC volume (mm³) | **262.94** | **174.99** |
| Difference PCCT − EID (mm³) | **87.95** |
| Mean of PCCT and EID (mm³) | **218.96** |

Old (canonical Gate 4 input) measurement for PT-129 reported on slide 9:
diff = 120.00 mm³, mean = 240.00 mm³.

## How to resume

**1. Re-authenticate AWS SSO (token TTL is short):**

```bash
aws sso login
aws sts get-caller-identity        # confirm role: AWSReservedSSO_SSMAccess_...
```

**2. Confirm venv + repo are still on the box** (just in case `/var/tmp` was cleared):

`C:/Users/EricaFreund/AppData/Local/Temp/ssm_setup3.json` already contains the
idempotent setup (creates venv if missing, clones repo if missing). Run it again:

```bash
export PYTHONUTF8=1; export PYTHONIOENCODING=utf-8
CID=$(aws ssm send-command \
  --cli-input-json "file://C:/Users/EricaFreund/AppData/Local/Temp/ssm_setup3.json" \
  --region us-east-1 \
  --query 'Command.CommandId' --output text)
echo "$CID"
# Then poll get-command-invocation until Status != InProgress|Pending.
# Read result by --output json > file then python json.load (avoids Windows charmap crash).
```

**3. Run on PT-129 first.** The wrapper script needs both workitem dirs:

```bash
WORK=/var/tmp/pcct_subseg
SCRIPTS=$WORK/PCCTvsEID/ComparePCCTandEIDWorkItems
OUT=$WORK/output/PT-129
mkdir -p $OUT
. $WORK/venv/bin/activate
bash $SCRIPTS/workitem_plaque_intersection_calculator.sh \
  --workitemA "/inst/zenith/AppData/Working Storage/wi_layne.cassidy/PT-129/wi-aff6fdcf" \
  --workitemB "/inst/zenith/AppData/Working Storage/wi_layne.cassidy/PT-129/wi-8622f0aa" \
  --output_dir $OUT
# Wrapper produces:
#   $OUT/final/right/common_vessel_sections_summary.json
#   $OUT/final/right/common_vessel_sections_volume_summary.txt
#   $OUT/final/left/common_vessel_sections_summary.json
#   $OUT/final/left/common_vessel_sections_volume_summary.txt
#   $OUT/inter/{right,left}/...
```

**4. Pull the JSON outputs back** (read-only `cat`, no scp needed):

```bash
aws ssm send-command --instance-ids i-097ce78e2d5f450ab \
  --document-name AWS-RunShellScript --region us-east-1 \
  --parameters 'commands=["cat /var/tmp/pcct_subseg/output/PT-129/final/right/common_vessel_sections_summary.json","echo ===LEFT===","cat /var/tmp/pcct_subseg/output/PT-129/final/left/common_vessel_sections_summary.json"]' \
  --query 'Command.CommandId' --output text
```

**5. Compare against the anchors above.**
- Match → proceed to the other 24 patients (loop the wrapper).
- Mismatch → inspect `common_vessel_sections_volume_summary.txt` for per-vessel
  breakdown; check that wrapper used the correct target side mapping for the
  workitem (LeftCoronary / RightCoronary).

## Key paths and IDs

| Thing | Where |
|---|---|
| PCCT framework root | `C:/Users/EricaFreund/OneDrive - Elucid Bioimaging Inc/PCCT/` |
| Current gate report (raw + length-normalized BA) | `PCCT/gate_results/gate_summary.txt` |
| Workitem summary CSVs (current Gate 3/4 input) | `PCCT/workitem_summaries/{PCCT,EID}/` |
| Target dir for sub-segment CSVs | `PCCT/workitem_summaries/subsegment/` (does not exist yet) |
| PT-129 verification deck | `PCCT/PCCT and CCTA intersection PT-129.pptx` |
| Dev box (us-east-1) | `i-097ce78e2d5f450ab` |
| SSM payload JSON files | `C:/Users/EricaFreund/AppData/Local/Temp/ssm_*.json` |
| SSM result JSON files | `C:/Users/EricaFreund/AppData/Local/Temp/*_result*.json` |

## Open SSM command IDs (for reference, may be expired)

| CID | Purpose | Status |
|---|---|---|
| 1fcd4f52-d260-4170-aacd-4c27eb4a328c | Scout PT-129 workitems | Done (Failed wrapper exit because system python lacked SimpleITK — but the find commands ran and gave us the paths) |
| ab2a5d60-09a4-4088-8a66-850e02ca7063 | Probe mounts | Success — confirmed /tmp noexec |
| 7cee7c3b-a7b3-435b-a243-0187fbcd28ec | First venv setup on /tmp | Failed — .so map error from noexec |
| (next) | Setup venv on /var/tmp | Was about to send `ssm_setup3.json` when paused |

## Constraints to keep in mind

- **Never write to `/inst/zenith/AppData/Working Storage`** — that's the shared workitem store.
- All scratch lives in `/var/tmp/pcct_subseg/` on the box; cleanup deferred until project complete.
- Canonical Gate 4 in `gate_summary.txt` stays unchanged; sub-segment Gate 3 + Gate 4 are added as **parallel sections**, not replacements.
- Windows AWS CLI: `aws ssm get-command-invocation` text output crashes on non-ASCII; always redirect `--output json` to a file and parse separately.
