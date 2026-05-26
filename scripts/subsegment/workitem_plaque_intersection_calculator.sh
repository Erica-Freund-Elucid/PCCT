#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  workitem_plaque_intersection_calculator.sh \
    --workitemA <dir> \
    --workitemB <dir> \
    --output_dir <output_dir>

This wrapper expects:
  <workitemA>/workitem.json
  <workitemB>/workitem.json

It runs common_vessel_sections_plaque_volume.py twice:
  1. right target
  2. left target

Pipeline call pattern:
  common_vessel_sections_plaque_volume.py \
    --study-a readings,lumen,wall,plaque_composition,wall_partition \
    --study-b readings,lumen,wall,plaque_composition,wall_partition \
    --target left|right \
    --out <output_dir>/final/<target> \
    --inter_out <output_dir>/inter \
    --volume-mode subvoxel \
    --min-common-length-mm 5.0 \
    --tolerance-mm 0.75
EOF
}

WORKITEM_A=""
WORKITEM_B=""
OUTPUT_DIR=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --workitemA)
      WORKITEM_A="${2:-}"
      shift 2
      ;;
    --workitemB)
      WORKITEM_B="${2:-}"
      shift 2
      ;;
    --output_dir)
      OUTPUT_DIR="${2:-}"
      shift 2
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "ERROR: Unknown argument: $1" >&2
      usage >&2
      exit 1
      ;;
  esac
done

if [[ -z "$WORKITEM_A" || -z "$WORKITEM_B" || -z "$OUTPUT_DIR" ]]; then
  echo "ERROR: --workitemA, --workitemB, and --output_dir are required." >&2
  usage >&2
  exit 1
fi

WORKITEM_A="$(realpath "$WORKITEM_A")"
WORKITEM_B="$(realpath "$WORKITEM_B")"
OUTPUT_DIR="$(realpath -m "$OUTPUT_DIR")"

WORKITEM_A_JSON="${WORKITEM_A}/workitem.json"
WORKITEM_B_JSON="${WORKITEM_B}/workitem.json"

if [[ ! -f "$WORKITEM_A_JSON" ]]; then
  echo "ERROR: Missing workitem JSON: $WORKITEM_A_JSON" >&2
  exit 1
fi

if [[ ! -f "$WORKITEM_B_JSON" ]]; then
  echo "ERROR: Missing workitem JSON: $WORKITEM_B_JSON" >&2
  exit 1
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_SCRIPT="${SCRIPT_DIR}/common_vessel_sections_plaque_volume.py"

if [[ ! -f "$PIPELINE_SCRIPT" ]]; then
  echo "ERROR: Pipeline script not found next to wrapper: $PIPELINE_SCRIPT" >&2
  exit 1
fi

extract_study_spec() {
  local workitem_dir="$1"
  local target="$2"

  python3 - "$workitem_dir" "$target" <<'PY'
import json
import sys
from pathlib import Path

workitem_dir = Path(sys.argv[1]).expanduser().resolve()
target = sys.argv[2].lower()
workitem_json = workitem_dir / "workitem.json"

def normalize_body_site(value):
    text = str(value or "").replace("_", " ").replace("-", " ").strip().lower()
    compact = "".join(ch for ch in text if ch.isalnum())

    if compact in {"leftcoronary", "leftcarotid", "left"}:
        return "left"
    if compact in {"rightcoronary", "rightcarotid", "right"}:
        return "right"

    if "left" in text and ("coronary" in text or "carotid" in text):
        return "left"
    if "right" in text and ("coronary" in text or "carotid" in text):
        return "right"

    return None

def get_nested(obj, keys, label):
    cur = obj
    for key in keys:
        if not isinstance(cur, dict) or key not in cur:
            raise KeyError(f"Missing {label}: {'.'.join(keys)}")
        cur = cur[key]
    if cur in (None, ""):
        raise ValueError(f"Empty {label}: {'.'.join(keys)}")
    return str(cur)

def abs_path(rel):
    p = Path(rel).expanduser()
    if p.is_absolute():
        return p.resolve()
    return (workitem_dir / p).resolve()

with workitem_json.open("r", encoding="utf-8") as f:
    payload = json.load(f)

target_defs = payload.get("targetDefinitions")
if not isinstance(target_defs, list):
    raise ValueError(f"{workitem_json} does not contain targetDefinitions array")

matched = []
for item in target_defs:
    if not isinstance(item, dict):
        continue
    side = normalize_body_site(item.get("bodySite"))
    if side == target:
        matched.append(item)

if not matched:
    raise ValueError(f"No targetDefinitions entry found for target={target} in {workitem_json}")

if len(matched) > 1:
    body_sites = [str(m.get("bodySite", "")) for m in matched]
    raise ValueError(
        f"Multiple targetDefinitions entries found for target={target} in {workitem_json}: {body_sites}"
    )

item = matched[0]

readings = abs_path(get_nested(item, ["readingsLocalFileName"], "readingsLocalFileName"))
lumen = abs_path(get_nested(item, ["regions", "lumenSegmentation"], "regions.lumenSegmentation"))
wall = abs_path(get_nested(item, ["regions", "wallSegmentation"], "regions.wallSegmentation"))
composition = abs_path(get_nested(item, ["probabilityMaps", "composition"], "probabilityMaps.composition"))
wall_partition = abs_path(get_nested(item, ["regions", "wallPartition"], "regions.wallPartition"))

paths = [readings, lumen, wall, composition, wall_partition]

missing = [str(p) for p in paths if not p.exists()]
if missing:
    raise FileNotFoundError(
        "Missing required files for target={} in {}:\n{}".format(
            target,
            workitem_json,
            "\n".join(missing),
        )
    )

print(",".join(str(p) for p in paths))
PY
}

run_for_target() {
  local target="$1"

  echo "============================================================"
  echo "Preparing ${target} target"
  echo "============================================================"

  local study_a
  local study_b

  study_a="$(extract_study_spec "$WORKITEM_A" "$target")"
  study_b="$(extract_study_spec "$WORKITEM_B" "$target")"

  local out_dir="${OUTPUT_DIR}/final/${target}"
  local inter_dir="${OUTPUT_DIR}/inter"

  mkdir -p "$out_dir" "$inter_dir"

  echo "Running pipeline for target=${target}"
  echo "Study A: ${study_a}"
  echo "Study B: ${study_b}"
  echo "Output:  ${out_dir}"
  echo "Inter:   ${inter_dir}"

  python3 "$PIPELINE_SCRIPT" \
    --study-a "$study_a" \
    --study-b "$study_b" \
    --target "$target" \
    --out "$out_dir" \
    --inter_out "$inter_dir" \
    --volume-mode subvoxel \
    --min-common-length-mm 5.0 \
    --tolerance-mm 0.75
}

run_for_target "right"
run_for_target "left"

echo "============================================================"
echo "Completed both targets"
echo "Final output:        ${OUTPUT_DIR}/final"
echo "Intermediate output: ${OUTPUT_DIR}/inter"
echo "============================================================"
