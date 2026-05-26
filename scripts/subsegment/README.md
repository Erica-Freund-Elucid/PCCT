# Common Vessel Plaque Volume Pipeline

This repository contains utilities to compare plaque, lumen, and wall volumes between two coronary/carotid studies over matched common vessel sections.

This compares resoluts in a PCCT and EID workitem of the same patient. It examines the two analyzed workitems and considers the plaque that lies in the intersection of the two workitem centerlines. Intersection is computed based on the vessel segment name and path length. It is granular enough to be sub-segment. ie half a D1 may be included in the intersection. The path length is computed based on the distance from the ostium.

This is equivalent to culling away the portion of one workitem that does not lie in the other workitem and then re-analyzing it. The resulting volumes in each workitem will be <= the original volumes.

The scripts are documented in this order:

1. [`run_calc_volume_pipeline.sh`](#1-runcalcvolumepipelinesh)
2. [`common_vessel_sections_plaque_volume.py`](#2-commonvesselsectionsplaquevolumepy)
3. [`multi_nrrd_reader_writer.py`](#3-multinrrdreaderwriterpy)

---

## Overview

The workflow compares two studies, commonly called **Study A** and **Study B**, for the same patient. Each study contains:

- readings JSON containing vessel centerline/cross-section information
- lumen SDF NRRD
- wall SDF NRRD
- plaque composition multi-NRRD
- wall partition multi-NRRD

The pipeline:

1. Reads workitem metadata or direct CLI inputs.
2. Extracts left and/or right target definitions.
3. Parses vessel centerlines from readings JSON.
4. Matches common vessel sections between the two studies.
5. Splits long vessels into discrete common subsections when matching is discontinuous.
6. Unpacks plaque and wall-partition multi-NRRD files.
7. Computes lumen, wall, CALC, LRNC, and NonCALCMATX volumes over common sections.
8. Writes JSON, TXT, VTP, SDF NRRD, and label-map QA outputs.

---

## Dependencies

Install the required Python packages in the environment used to run the scripts:

```bash
pip install numpy scipy SimpleITK pynrrd vtk
```

Notes:

- `SimpleITK` is required for reading/writing NRRD, NIfTI, and signed distance maps.
- `pynrrd` is required by `multi_nrrd_reader_writer.py`.
- `scipy` is used for faster nearest-neighbor distance queries.
- `vtk` is required when using `--volume-mode subvoxel`.

The wrapper script also requires:

```bash
bash
python3
realpath
```

---

# 1. `run_calc_volume_pipeline.sh`

## Purpose

`run_calc_volume_pipeline.sh` is the recommended top-level entry point when each study is represented by a workitem directory containing a `workitem.json`.

The wrapper:

1. Accepts two workitem directories.
2. Reads each `workitem.json`.
3. Extracts file paths for both the **right** and **left** target definitions.
4. Builds the compact `--study-a` and `--study-b` arguments required by `common_vessel_sections_plaque_volume.py`.
5. Runs the Python pipeline twice:
   - first for `right`
   - then for `left`

The script expects `common_vessel_sections_plaque_volume.py` to be in the same directory as the shell script.

## Usage

```bash
./run_calc_volume_pipeline.sh \
  --workitemA <workitem_a_dir> \
  --workitemB <workitem_b_dir> \
  --output_dir <output_dir>
```

## Parameters

| Parameter | Required | Description |
|---|---:|---|
| `--workitemA <dir>` | Yes | Directory containing Study A `workitem.json`. |
| `--workitemB <dir>` | Yes | Directory containing Study B `workitem.json`. |
| `--output_dir <dir>` | Yes | Root output directory. The wrapper creates `final/right`, `final/left`, and `inter` beneath this directory. |
| `-h`, `--help` | No | Prints usage information. |

## Expected `workitem.json` structure

Each workitem directory must contain:

```text
<workitem_dir>/workitem.json
```

The JSON must contain a `targetDefinitions` array. For each target definition, the wrapper reads:

```text
bodySite
readingsLocalFileName
regions.lumenSegmentation
regions.wallSegmentation
probabilityMaps.composition
regions.wallPartition
```

The wrapper maps body sites to target sides as follows:

| `bodySite` examples | Pipeline target |
|---|---|
| `LeftCoronary`, `Left Carotid`, `left` | `left` |
| `RightCoronary`, `RightCarotid`, `right` | `right` |

All file paths in the workitem JSON are treated as relative to the workitem directory unless already absolute.

## Internal Python calls

For the right target, the wrapper calls:

```bash
python3 common_vessel_sections_plaque_volume.py \
  --study-a <study_a_right_spec> \
  --study-b <study_b_right_spec> \
  --target right \
  --out <output_dir>/final/right \
  --inter_out <output_dir>/inter \
  --volume-mode subvoxel \
  --min-common-length-mm 5.0 \
  --tolerance-mm 0.75
```

For the left target, the wrapper calls:

```bash
python3 common_vessel_sections_plaque_volume.py \
  --study-a <study_a_left_spec> \
  --study-b <study_b_left_spec> \
  --target left \
  --out <output_dir>/final/left \
  --inter_out <output_dir>/inter \
  --volume-mode subvoxel \
  --min-common-length-mm 5.0 \
  --tolerance-mm 0.75
```

## Output layout

Given:

```bash
--output_dir out_dir
```

the wrapper writes:

```text
out_dir/
  final/
    right/
      common_vessel_sections_summary.json
      common_vessel_sections_volume_summary.txt
      common_vessel_sections_geometry_summary.txt
    left/
      common_vessel_sections_summary.json
      common_vessel_sections_volume_summary.txt
      common_vessel_sections_geometry_summary.txt
  inter/
    right/
      ...
    left/
      ...
```

## Sample command

```bash
./run_calc_volume_pipeline.sh \
  --workitemA /data/patient001/pcct \
  --workitemB /data/patient001/ccta \
  --output_dir /data/patient001/calc_volume_output
```


---

# 2. `common_vessel_sections_plaque_volume.py`

## Purpose

`common_vessel_sections_plaque_volume.py` is the main pipeline script. It compares two studies for either the left or right coronary/carotid target and computes common-section volumes for:

- lumen
- wall
- CALC
- LRNC
- NonCALCMATX

It also writes intermediate QA artifacts for centerlines and restricted component SDFs.

## High-level processing logic

### Unit 1: common vessel section matching

1. Parses recursive readings JSON.
2. Extracts vessel cross-sections including:
   - `position`
   - `xaxis`
   - `yaxis`
   - `path_distance`
   - `segment_distance`
   - `vessel_distance`
   - vessel/segment name
3. Preserves proximal-to-distal JSON order.
4. Normalizes vessel names using aliases such as:
   - `LAD` -> `LeftAnteriorDescending`
   - `LCX` -> `Circumflex`
   - `LM` -> `MainStem`
5. Finds vessels present in both studies.
6. Matches cross-section samples using distance tolerance.
7. Splits long vessels into multiple discrete common subsections when common matching is discontinuous.
8. Writes diagnostic VTP centerlines.

### Unit 2: SDF volume processing

1. Unpacks plaque composition multi-NRRD using `multi_nrrd_reader_writer.py`.
2. Unpacks wall partition multi-NRRD using `multi_nrrd_reader_writer.py`.
3. For each common vessel subsection:
   - selects the vessel-specific wall partition SDF
   - clips using proximal and distal cross-section planes
   - applies a centerline safety radius
   - computes component volumes inside the common section
4. Writes:
   - summary JSON
   - volume comparison TXT report
   - common-section geometry TXT report
   - per-section and total SDF NRRD QA outputs
   - combined per-study NIfTI label maps

## Usage

```bash
python3 common_vessel_sections_plaque_volume.py \
  --study-a <readings,lumen,wall,plaque_multi,wall_partition_multi> \
  --study-b <readings,lumen,wall,plaque_multi,wall_partition_multi> \
  --target <left|right> \
  --out <output_dir> \
  --inter_out <intermediate_output_dir> \
  [options]
```

## Required parameters

| Parameter | Required | Description |
|---|---:|---|
| `--study-a` | Yes | Comma-separated Study A inputs in fixed order: `readings.json,lumen.nrrd,wall.nrrd,plaque.multi.nrrd,wallPartition.multi.nrrd`. |
| `--study-b` | Yes | Comma-separated Study B inputs in fixed order: `readings.json,lumen.nrrd,wall.nrrd,plaque.multi.nrrd,wallPartition.multi.nrrd`. |
| `--target {left,right}` | Yes | Target vascular tree side. |
| `--out <dir>` | Yes | Final output directory. |
| `--inter_out <dir>` | Yes | Intermediate/debug output directory. |

## Optional parameters

| Parameter | Default | Description |
|---|---:|---|
| `--volume-mode {voxel,subvoxel}` | `voxel` | Volume calculation mode. `voxel` uses voxel count times voxel volume. `subvoxel` extracts the zero-level surface with VTK and computes mesh volume. |
| `--sdf-threshold-mm <float>` | `0.0` | SDF inside threshold. Default means voxels with `SDF <= 0` are inside. |
| `--centerline-radius-mm <float>` | `25.0` | Safety radius around common centerline samples. |
| `--cap-tolerance-mm <float>` | `0.5` | Tolerance for proximal/distal cross-section cap planes. |
| `--tolerance-mm <float>` | `0.75` | Centerline distance tolerance for matching common cross-section samples. |
| `--min-common-length-mm <float>` | `5.0` | Minimum accepted common subsection length. |
| `--max-common-gap-mm <float>` | Derived automatically | Maximum allowed gap within a common subsection. If omitted, it is derived from cross-section sampling and tolerance. |

## Study argument format

Each study argument must contain exactly five comma-separated paths:

```text
readings.json,lumen.nrrd,wall.nrrd,plaque.multi.nrrd,wallPartition.multi.nrrd
```

Example:

```text
/data/pcct/readings.json,/data/pcct/lumen.nrrd,/data/pcct/wall.nrrd,/data/pcct/composition.multi.nrrd,/data/pcct/wallPartition.multi.nrrd
```

## Volume modes

### `voxel`

Voxel mode computes:

```text
volume_mm3 = voxel_count * spacing_x * spacing_y * spacing_z
```

### `subvoxel`

Sub-voxel mode computes the zero-level-set surface using VTK:

```text
vtkFlyingEdges3D -> vtkTriangleFilter -> vtkCleanPolyData -> vtkMassProperties
```

For plaque or other component volumes inside the common vessel section, level-set intersection is used:

```text
component_common_sdf = max(component_sdf, common_section_sdf)
```

For total volume over multiple common subsections, level-set union is used:

```text
total_component_sdf = min(section_1_sdf, section_2_sdf, ...)
```

## Main outputs

### Final outputs

Written under `--out`:

```text
common_vessel_sections_summary.json
common_vessel_sections_volume_summary.txt
common_vessel_sections_geometry_summary.txt
```

`common_vessel_sections_summary.json` contains input paths, target, volume mode, matching parameters, generated file paths, matched vessel summaries, and volume results.

`common_vessel_sections_volume_summary.txt` contains tabular volume results for:

- Wall
- Lumen
- Plaque - CALC
- Plaque - LRNC
- Plaque - NonCALCMATX

It includes vessel totals, subsection values, Study A/Study B values, differences, and percent deltas for total rows.

`common_vessel_sections_geometry_summary.txt` contains vessel aggregate common length, subsection length, and Study A/B start/end XYZ coordinates.

### Intermediate outputs

Written under `--inter_out/<target>`:

```text
common_vessel_sections.vtp
study_a_complete_vessels.vtp
study_b_complete_vessels.vtp
```

Additional outputs are written under:

```text
--inter_out/<target>/study_a/
--inter_out/<target>/study_b/
```

These include:

```text
plaque/
wallPartition/
common_sections/
```

## Sample direct command

```bash
python3 common_vessel_sections_plaque_volume.py \
  --study-a /data/pcct/readings.json,/data/pcct/lumen.nrrd,/data/pcct/wall.nrrd,/data/pcct/composition.multi.nrrd,/data/pcct/wallPartition.multi.nrrd \
  --study-b /data/ccta/readings.json,/data/ccta/lumen.nrrd,/data/ccta/wall.nrrd,/data/ccta/composition.multi.nrrd,/data/ccta/wallPartition.multi.nrrd \
  --target left \
  --out /data/output/final/left \
  --inter_out /data/output/inter \
  --volume-mode subvoxel \
  --min-common-length-mm 5.0 \
  --tolerance-mm 0.75
```


---

# 3. `multi_nrrd_reader_writer.py`

## Purpose

`multi_nrrd_reader_writer.py` reads and writes 4D `.multi.nrrd` files.

It is used by the main pipeline to unpack:

- plaque composition multi-NRRD
- wall partition multi-NRRD

into individual 3D SDF NRRD files.

The script supports both:

1. region-based metadata, such as:
   - `slice0_region`
   - `slice1_region`
2. site-based metadata, such as:
   - `slice0_site`
   - `slice1_site`

If site metadata is present, individual output files are named by site. Otherwise, they are named by region.

## Modes

The script has two modes:

```text
read
write
```

---

## Read mode

### Purpose

Read mode splits a 4D `.multi.nrrd` into individual 3D `.nrrd` files.

### Usage

```bash
python3 multi_nrrd_reader_writer.py read \
  --input <input.multi.nrrd> \
  --output-dir <output_dir>
```

### Parameters

| Parameter | Required | Description |
|---|---:|---|
| `--input <file>` | Yes | Input 4D `.multi.nrrd`. |
| `--output-dir <dir>` | Yes | Directory where individual 3D `.nrrd` files are written. |

### Output naming

If the input header has site metadata:

```text
slice0_site
slice1_site
...
```

then files are written by site name:

```text
<site>.nrrd
```

Otherwise, files are written by region name:

```text
CALC.nrrd
LRNC.nrrd
NonCALCMATX.nrrd
```

### Example: unpack plaque composition

```bash
python3 multi_nrrd_reader_writer.py read \
  --input /data/pcct/composition.multi.nrrd \
  --output-dir /data/output/inter/left/study_a/plaque
```

Example outputs:

```text
/data/output/inter/left/study_a/plaque/CALC.nrrd
/data/output/inter/left/study_a/plaque/LRNC.nrrd
/data/output/inter/left/study_a/plaque/NonCALCMATX.nrrd
```

### Example: unpack wall partition

```bash
python3 multi_nrrd_reader_writer.py read \
  --input /data/pcct/wallPartition.multi.nrrd \
  --output-dir /data/output/inter/left/study_a/wallPartition
```

Example outputs:

```text
/data/output/inter/left/study_a/wallPartition/LeftAnteriorDescending.nrrd
/data/output/inter/left/study_a/wallPartition/Circumflex.nrrd
/data/output/inter/left/study_a/wallPartition/MainStem.nrrd
```

---

## Write mode

### Purpose

Write mode combines individual 3D NRRDs into a 4D `.multi.nrrd`.

### Usage

```bash
python3 multi_nrrd_reader_writer.py write \
  --input-dir <input_dir> \
  --output-dir <output_dir> \
  --target-path <target_path> \
  [options]
```

### Parameters

| Parameter | Required | Description |
|---|---:|---|
| `--input-dir <dir>` | Yes | Directory containing individual 3D NRRD files. |
| `--output-dir <dir>` | Yes | Directory where the 4D `.multi.nrrd` is written. |
| `--target-path <target_path>` | Yes | Target path metadata, for example `LeftCoronary`. |
| `--output-name <name>` | No | Output multi-NRRD filename. Default is `composition.multi.nrrd`. |
| `--sites <site1> <site2> ...` | No | Optional site names. When provided, input files are read as `<site>.nrrd`, and `slice{i}_site` metadata is written. |
| `--slice-metadata-json <file>` | No | Optional JSON describing slice metadata. Supports `region`, `site`, `pixelInterpretation`, and `targetPath`. |

### Region-based write example

```bash
python3 multi_nrrd_reader_writer.py write \
  --input-dir /data/split_plaque \
  --output-dir /data/repacked \
  --target-path LeftCoronary \
  --output-name composition.multi.nrrd
```

Expected input files:

```text
/data/split_plaque/CALC.nrrd
/data/split_plaque/LRNC.nrrd
/data/split_plaque/IPH.nrrd
/data/split_plaque/PVAT.nrrd
/data/split_plaque/MATX.nrrd
/data/split_plaque/FIBL.nrrd
/data/split_plaque/NonCALCMATX.nrrd
```

### Site-aware write example

```bash
python3 multi_nrrd_reader_writer.py write \
  --input-dir /data/split_wall_partition \
  --output-dir /data/repacked \
  --target-path LeftCoronary \
  --output-name wallPartition.multi.nrrd \
  --sites LeftAnteriorDescending Circumflex LeftMarginal1
```

Expected input files:

```text
/data/split_wall_partition/LeftAnteriorDescending.nrrd
/data/split_wall_partition/Circumflex.nrrd
/data/split_wall_partition/LeftMarginal1.nrrd
```

The output multi-NRRD header will include:

```text
slice0_site:=LeftAnteriorDescending
slice1_site:=Circumflex
slice2_site:=LeftMarginal1
```

### Slice metadata JSON example

```json
{
  "slices": [
    {
      "region": "CALC",
      "site": "ProximalLAD",
      "pixelInterpretation": "SignedDistanceFunction",
      "targetPath": "LeftCoronary"
    },
    {
      "region": "LRNC",
      "site": "MidLAD",
      "pixelInterpretation": "SignedDistanceFunction",
      "targetPath": "LeftCoronary"
    }
  ]
}
```

Run:

```bash
python3 multi_nrrd_reader_writer.py write \
  --input-dir /data/split_nrrds \
  --output-dir /data/repacked \
  --target-path LeftCoronary \
  --slice-metadata-json /data/slice_metadata.json
```

---

## Typical full workflow

### Option 1: use the wrapper

```bash
./run_calc_volume_pipeline.sh \
  --workitemA /data/patient001/pcct \
  --workitemB /data/patient001/ccta \
  --output_dir /data/patient001/output
```

### Option 2: call the main pipeline directly

```bash
python3 common_vessel_sections_plaque_volume.py \
  --study-a /data/pcct/readings.json,/data/pcct/lumen.nrrd,/data/pcct/wall.nrrd,/data/pcct/composition.multi.nrrd,/data/pcct/wallPartition.multi.nrrd \
  --study-b /data/ccta/readings.json,/data/ccta/lumen.nrrd,/data/ccta/wall.nrrd,/data/ccta/composition.multi.nrrd,/data/ccta/wallPartition.multi.nrrd \
  --target left \
  --out /data/output/final/left \
  --inter_out /data/output/inter \
  --volume-mode subvoxel \
  --min-common-length-mm 5.0 \
  --tolerance-mm 0.75
```

---

---

## Output directory - files and sub-directory structure

out_dir/
├── final/
│   ├── left/
│   │   ├── common_vessel_sections_geometry_summary.txt
│   │   ├── common_vessel_sections_summary.json
│   │   └── common_vessel_sections_volume_summary.txt
│   │
│   └── right/
│       ├── common_vessel_sections_geometry_summary.txt
│       ├── common_vessel_sections_summary.json
│       └── common_vessel_sections_volume_summary.txt
│
└── inter/
    ├── left/
    │   ├── common_vessel_sections.vtp
    │   ├── study_a_complete_vessels.vtp
    │   ├── study_b_complete_vessels.vtp
    │   │
    │   ├── study_a/
    │   │   ├── plaque/
    │   │   │   ├── \<plaque_component_1\>.nrrd
    │   │   │   ├── \<plaque_component_2\>.nrrd
    │   │   │   └── ...
    │   │   │
    │   │   ├── wallPartition/
    │   │   │   ├── \<vessel_wall_partition_1\>.nrrd
    │   │   │   ├── \<vessel_wall_partition_2\>.nrrd
    │   │   │   └── ...
    │   │   │
    │   │   └── common_sections/
    │   │       ├── study_a_combined_common_sections_label.nii.gz
    │   │       │
    │   │       ├── \<common_vessel_1\>/
    │   │       │   ├── section_001/
    │   │       │   │   ├── \<plaque_component_1\>_common.sdf.nrrd
    │   │       │   │   ├── \<plaque_component_2\>_common.sdf.nrrd
    │   │       │   │   └── ...
    │   │       │   │
    │   │       │   ├── section_002/
    │   │       │   │   ├── \<plaque_component_1>_common.sdf.nrrd
    │   │       │   │   └── ...
    │   │       │   │
    │   │       │   ├── ...
    │   │       │   │
    │   │       │   └── total/
    │   │       │       ├── \<plaque_component_1\>_common_total.sdf.nrrd
    │   │       │       └── ...
    │   │       │
    │   │       ├── \<common_vessel_2\>/
    │   │       │   ├── ...
    │   │       │
    │   │       └── ...
    │   │
    │   └── study_b/
    │       ├── plaque/
    │       │   ├── \<plaque_component_1\>.nrrd
    │       │   ├── \<plaque_component_2\>.nrrd
    │       │   └── ...
    │       │
    │       ├── wallPartition/
    │       │   ├── \<vessel_wall_partition_1\>.nrrd
    │       │   ├── \<vessel_wall_partition_2\>.nrrd
    │       │   └── ...
    │       │
    │       └── common_sections/
    │           ├── study_b_combined_common_sections_label.nii.gz
    │           │
    │           ├── \<common_vessel_1\>/
    │           │   ├── section_001/
    │           │   │   ├── CALC_common.sdf.nrrd
    │           │   │   ├── LRNC_common.sdf.nrrd
    │           │   │   ├── NonCALCMATX_common.sdf.nrrd
    │           │   │   ├── lumen_common.sdf.nrrd
    │           │   │   └── wall_common.sdf.nrrd
    │           │   ├── section_002/
    │           │   │   ├── CALC_common.sdf.nrrd
    │           │   │   ├── LRNC_common.sdf.nrrd
    │           │   │   ├── NonCALCMATX_common.sdf.nrrd
    │           │   │   ├── lumen_common.sdf.nrrd
    │           │   │   └── wall_common.sdf.nrrd
    │           │   ├── ...
    │           │   └── total/
    │           │       ├── CALC_common_total.sdf.nrrd
    │           │       ├── LRNC_common_total.sdf.nrrd
    │           │       ├── NonCALCMATX_common_total.sdf.nrrd
    │           │       ├── lumen_common_total.sdf.nrrd
    │           │       └── wall_common_total.sdf.nrrd
    │           │
    │           ├── \<common_vessel_2\>/
    │           │   ├── section_001/
    │           │   │   ├── CALC_common.sdf.nrrd
    │           │   │   ├── LRNC_common.sdf.nrrd
    │           │   │   ├── NonCALCMATX_common.sdf.nrrd
    │           │   │   ├── lumen_common.sdf.nrrd
    │           │   │   └── wall_common.sdf.nrrd
    │           │   └── total/
    │           │       ├── CALC_common_total.sdf.nrrd
    │           │       ├── LRNC_common_total.sdf.nrrd
    │           │       ├── NonCALCMATX_common_total.sdf.nrrd
    │           │       ├── lumen_common_total.sdf.nrrd
    │           │       └── wall_common_total.sdf.nrrd
    │           │
    │           └── ...
    │
    └── right/
        └── ... Similar structure as above
---

## Troubleshooting

### `Pipeline script not found next to wrapper`

Make sure these files are in the same directory:

```text
run_calc_volume_pipeline.sh
common_vessel_sections_plaque_volume.py
multi_nrrd_reader_writer.py
```

### `No targetDefinitions entry found`

Check that `workitem.json` contains a `targetDefinitions` entry with `bodySite` matching the requested side:

```text
LeftCoronary, LeftCarotid, RightCoronary, RightCarotid
```

### `Missing required files`

The wrapper validates all file paths derived from `workitem.json`. Confirm the referenced files exist relative to the workitem directory, or use absolute paths in the JSON.

### `VTK is required for --volume-mode subvoxel`

Install VTK:

```bash
pip install vtk
```

### `SimpleITK is required`

Install SimpleITK:

```bash
pip install SimpleITK
```

### `pynrrd` import error

Install pynrrd:

```bash
pip install pynrrd
```
