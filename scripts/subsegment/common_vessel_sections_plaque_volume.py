#!/usr/bin/env python3
"""
Pipeline for comparing two coronary/carotid readings studies.

Unit 1:
  - Parse recursive readings JSON centerline/cross-sections.
  - Match common vessels between two studies.
  - Split common vessel paths into discrete subsections.
  - Write diagnostic VTPs.

Unit 2:
  - Inputs use compact fixed-order study specs:
      --study-a readings.json,lumen.nrrd,wall.nrrd,plaque.multi.nrrd,wallPartition.multi.nrrd
      --study-b readings.json,lumen.nrrd,wall.nrrd,plaque.multi.nrrd,wallPartition.multi.nrrd
  - Unpack plaque.multi.nrrd and wallPartition.multi.nrrd using multi_nrrd_reader_writer.py
    from the same directory as this script.
  - For each common vessel subsection, use the corresponding vessel wallPartition SDF,
    plus proximal/distal cross-section cap planes from xaxis/yaxis, to compute lumen,
    wall, CALC, LRNC, and NonCALCMATX volumes.
  - Write restricted component label maps (*.nii.gz) and SDFs (*.sdf.nrrd) for QA.
"""

from __future__ import annotations

import argparse
import json
import math
import re
import subprocess
import sys
import xml.etree.ElementTree as ET
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np

try:
    import SimpleITK as sitk
except ImportError:
    sitk = None

try:
    from scipy.spatial import cKDTree
except ImportError:
    cKDTree = None

try:
    import vtk
    from vtk.util import numpy_support as vtk_numpy_support
except ImportError:
    vtk = None
    vtk_numpy_support = None



def expand_path(value: str | Path) -> Path:
    """
    Expand shell-style user paths such as ~/data/file.nrrd.

    This is required because subprocess.run([...]) does not perform shell
    expansion when shell=False, which is the safer/default mode.
    """
    return Path(value).expanduser().resolve()

# -----------------------------
# Data model
# -----------------------------

@dataclass
class CrossSection:
    vessel_name: str
    segment_name: str
    position: Tuple[float, float, float]
    xaxis: Optional[Tuple[float, float, float]]
    yaxis: Optional[Tuple[float, float, float]]
    path_distance: Optional[float]
    segment_distance: Optional[float]
    vessel_distance: float
    original_index: int
    tree_path: str


@dataclass
class CommonSubSectionSummary:
    subsection_id: int
    study_a_points_common: int
    study_b_points_common: int
    common_start_mm: float
    common_end_mm: float
    common_length_mm: float


@dataclass
class VesselMatchSummary:
    vessel_name: str
    vessel_id: int
    study_a_points_total: int
    study_b_points_total: int
    study_a_points_common_total: int
    study_b_points_common_total: int
    total_common_length_mm: float
    subsections: List[CommonSubSectionSummary]


@dataclass
class StudyInputs:
    readings: Path
    lumen: Path
    wall: Path
    plaque_multi: Path
    wall_partition_multi: Path


# -----------------------------
# Readings JSON parsing
# -----------------------------

CANONICAL_ALIASES = {
    "mainstem": "MainStem",
    "leftmain": "MainStem",
    "lm": "MainStem",
    "lad": "LeftAnteriorDescending",
    "leftanteriordescending": "LeftAnteriorDescending",
    "left_anterior_descending": "LeftAnteriorDescending",
    "lcx": "Circumflex",
    "cx": "Circumflex",
    "circumflex": "Circumflex",
    "diagonal1": "Diagonal1",
    "d1": "Diagonal1",
    "firstdiagonal": "Diagonal1",
    "om1": "OM1",
    "obtusemarginal1": "OM1",
    "firstobtusemarginal": "OM1",
}


def normalize_name(name: Optional[str]) -> str:
    if not name:
        return ""
    cleaned = re.sub(r"[^A-Za-z0-9_]+", "", str(name)).strip()
    return CANONICAL_ALIASES.get(cleaned.lower(), cleaned)


def infer_name_from_cross_section(cs: Dict[str, Any]) -> str:
    for key in ("vessel_name", "vesselName", "segment_name", "segmentName", "name"):
        value = normalize_name(cs.get(key))
        if value:
            return value
    return ""


def as_vec3(value: Any) -> Optional[Tuple[float, float, float]]:
    if not isinstance(value, Sequence) or len(value) != 3:
        return None
    try:
        return (float(value[0]), float(value[1]), float(value[2]))
    except (TypeError, ValueError):
        return None


def numeric_or_none(x: Any) -> Optional[float]:
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def distance(a: Tuple[float, float, float], b: Tuple[float, float, float]) -> float:
    return math.sqrt((a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2 + (a[2] - b[2]) ** 2)


def cumulative_distances(points: Sequence[Tuple[float, float, float]]) -> List[float]:
    if not points:
        return []
    out = [0.0]
    for i in range(1, len(points)):
        out.append(out[-1] + distance(points[i - 1], points[i]))
    return out


def extract_root_segment(doc: Dict[str, Any]) -> Dict[str, Any]:
    if isinstance(doc.get("root_segment"), dict):
        return doc["root_segment"]
    if isinstance(doc.get("initialization_points"), dict):
        return doc["initialization_points"]
    raise ValueError("Could not find 'root_segment' or 'initialization_points' in readings JSON.")


def collect_segments_recursive(
    segment: Dict[str, Any],
    inherited_name: str = "",
    tree_path: str = "root",
) -> List[Tuple[str, str, List[Dict[str, Any]], str]]:
    cross_sections = segment.get("cross_sections") or []
    segment_name = normalize_name(segment.get("segment_name") or segment.get("name"))
    local_names = [infer_name_from_cross_section(cs) for cs in cross_sections]
    explicit_names = [n for n in local_names if n]

    if explicit_names:
        vessel_name = max(set(explicit_names), key=explicit_names.count)
    else:
        vessel_name = segment_name or inherited_name or "Unknown"

    vessel_name = normalize_name(vessel_name) or "Unknown"
    segment_name = segment_name or vessel_name

    chunks = []
    if cross_sections:
        chunks.append((vessel_name, segment_name, cross_sections, tree_path))

    for idx, child in enumerate(segment.get("distal_segments") or []):
        if isinstance(child, dict):
            chunks.extend(
                collect_segments_recursive(
                    child,
                    inherited_name=vessel_name,
                    tree_path=f"{tree_path}/distal_segments[{idx}]",
                )
            )
    return chunks


def choose_distance_source(raw_cross_sections: Sequence[Dict[str, Any]]) -> str:
    for key in ("vessel_distance", "path_distance", "segment_distance"):
        vals = [numeric_or_none(cs.get(key)) for cs in raw_cross_sections]
        vals = [v for v in vals if v is not None]
        if len(vals) < 2:
            continue
        if max(vals) - min(vals) <= 1e-6:
            continue
        if all(vals[i] >= vals[i - 1] for i in range(1, len(vals))):
            return key
    return "geometric"


def parse_readings_json(path: Path) -> Dict[str, List[CrossSection]]:
    with path.open("r", encoding="utf-8") as f:
        doc = json.load(f)

    root = extract_root_segment(doc)
    chunks = collect_segments_recursive(root)

    vessel_map: Dict[str, List[CrossSection]] = {}
    global_index = 0
    vessel_last_distance: Dict[str, float] = {}
    vessel_last_position: Dict[str, Tuple[float, float, float]] = {}

    for vessel_name, segment_name, raw_cross_sections, tree_path in chunks:
        positions = [as_vec3(cs.get("position")) for cs in raw_cross_sections]
        valid_positions = [p for p in positions if p is not None]
        geom_dist = cumulative_distances(valid_positions)
        source = choose_distance_source(raw_cross_sections)

        if source == "geometric":
            raw_distances = geom_dist
        else:
            raw_distances = []
            for cs in raw_cross_sections:
                val = numeric_or_none(cs.get(source))
                raw_distances.append(0.0 if val is None else float(val))
            if raw_distances:
                first = raw_distances[0]
                raw_distances = [d - first for d in raw_distances]

        if vessel_name in vessel_last_distance:
            offset = vessel_last_distance[vessel_name]
            if valid_positions and vessel_name in vessel_last_position:
                offset += distance(vessel_last_position[vessel_name], valid_positions[0])
        else:
            offset = 0.0

        chunk_sections: List[CrossSection] = []
        valid_cursor = 0
        for local_idx, cs in enumerate(raw_cross_sections):
            pos = as_vec3(cs.get("position"))
            if pos is None:
                continue

            if valid_cursor < len(raw_distances):
                local_d = raw_distances[valid_cursor]
            elif valid_cursor < len(geom_dist):
                local_d = geom_dist[valid_cursor]
            else:
                local_d = 0.0

            explicit = infer_name_from_cross_section(cs)
            chunk_sections.append(
                CrossSection(
                    vessel_name=vessel_name,
                    segment_name=explicit or segment_name,
                    position=pos,
                    xaxis=as_vec3(cs.get("xaxis")),
                    yaxis=as_vec3(cs.get("yaxis")),
                    path_distance=numeric_or_none(cs.get("path_distance")),
                    segment_distance=numeric_or_none(cs.get("segment_distance")),
                    vessel_distance=float(offset + local_d),
                    original_index=global_index,
                    tree_path=f"{tree_path}/cross_sections[{local_idx}]",
                )
            )
            global_index += 1
            valid_cursor += 1

        if chunk_sections:
            # Preserve proximal-to-distal JSON order. Do not sort by distance.
            vessel_map.setdefault(vessel_name, []).extend(chunk_sections)
            vessel_last_distance[vessel_name] = chunk_sections[-1].vessel_distance
            vessel_last_position[vessel_name] = chunk_sections[-1].position

    return vessel_map


# -----------------------------
# Common section matching
# -----------------------------

def build_vessel_ids(vessels_a: Dict[str, List[CrossSection]], vessels_b: Dict[str, List[CrossSection]]) -> Dict[str, int]:
    return {name: idx for idx, name in enumerate(sorted(set(vessels_a) | set(vessels_b)), start=1)}


def matching_vessel_names(vessels_a: Dict[str, List[CrossSection]], vessels_b: Dict[str, List[CrossSection]]) -> List[str]:
    return sorted(set(vessels_a) & set(vessels_b))


def nearest_distance_match_index(d: float, other: Sequence[CrossSection], tolerance_mm: float) -> Optional[int]:
    if not other:
        return None
    distances = [abs(d - s.vessel_distance) for s in other]
    idx = int(np.argmin(distances))
    return idx if distances[idx] <= tolerance_mm else None


def split_common_subsections(
    a: Sequence[CrossSection],
    b: Sequence[CrossSection],
    tolerance_mm: float,
    min_common_length_mm: float,
    max_common_gap_mm: Optional[float],
) -> List[Tuple[List[CrossSection], List[CrossSection]]]:
    if len(a) < 2 or len(b) < 2:
        return []

    if max_common_gap_mm is None:
        a_gaps = [a[i].vessel_distance - a[i - 1].vessel_distance for i in range(1, len(a)) if a[i].vessel_distance >= a[i - 1].vessel_distance]
        typical_gap = float(np.median(a_gaps)) if a_gaps else tolerance_mm
        max_common_gap_mm = max(2.0 * tolerance_mm, 3.0 * typical_gap)

    matched_pairs: List[Tuple[int, int]] = []
    for ai, section in enumerate(a):
        bi = nearest_distance_match_index(section.vessel_distance, b, tolerance_mm)
        if bi is not None:
            matched_pairs.append((ai, bi))

    if not matched_pairs:
        return []

    runs: List[List[Tuple[int, int]]] = []
    current = [matched_pairs[0]]
    for ai, bi in matched_pairs[1:]:
        prev_ai, prev_bi = current[-1]
        a_gap = a[ai].vessel_distance - a[prev_ai].vessel_distance
        b_gap = b[bi].vessel_distance - b[prev_bi].vessel_distance
        if ai > prev_ai and bi >= prev_bi and a_gap <= max_common_gap_mm and b_gap <= max_common_gap_mm:
            current.append((ai, bi))
        else:
            runs.append(current)
            current = [(ai, bi)]
    runs.append(current)

    subsections: List[Tuple[List[CrossSection], List[CrossSection]]] = []
    for run in runs:
        if len(run) < 2:
            continue
        a_idx = [p[0] for p in run]
        b_idx = [p[1] for p in run]
        a_sections = list(a[min(a_idx): max(a_idx) + 1])
        b_sections = list(b[min(b_idx): max(b_idx) + 1])
        if len(a_sections) < 2 or len(b_sections) < 2:
            continue
        length = max(0.0, min(
            a_sections[-1].vessel_distance - a_sections[0].vessel_distance,
            b_sections[-1].vessel_distance - b_sections[0].vessel_distance,
        ))
        if length + 1e-6 >= min_common_length_mm:
            subsections.append((a_sections, b_sections))
    return subsections


def compute_common_sections(
    vessels_a: Dict[str, List[CrossSection]],
    vessels_b: Dict[str, List[CrossSection]],
    tolerance_mm: float,
    vessel_ids: Dict[str, int],
    min_common_length_mm: float,
    max_common_gap_mm: Optional[float],
) -> Tuple[Dict[str, List[Tuple[List[CrossSection], List[CrossSection]]]], List[VesselMatchSummary]]:
    common: Dict[str, List[Tuple[List[CrossSection], List[CrossSection]]]] = {}
    summaries: List[VesselMatchSummary] = []

    for vessel_name in matching_vessel_names(vessels_a, vessels_b):
        a = vessels_a[vessel_name]
        b = vessels_b[vessel_name]
        subsections = split_common_subsections(a, b, tolerance_mm, min_common_length_mm, max_common_gap_mm)
        if not subsections:
            continue

        common[vessel_name] = subsections
        sub_summaries: List[CommonSubSectionSummary] = []
        total_a = 0
        total_b = 0
        total_len = 0.0
        for i, (a_sec, b_sec) in enumerate(subsections, start=1):
            start_mm = max(a_sec[0].vessel_distance, b_sec[0].vessel_distance)
            end_mm = min(a_sec[-1].vessel_distance, b_sec[-1].vessel_distance)
            length_mm = max(0.0, end_mm - start_mm)
            total_a += len(a_sec)
            total_b += len(b_sec)
            total_len += length_mm
            sub_summaries.append(CommonSubSectionSummary(i, len(a_sec), len(b_sec), start_mm, end_mm, length_mm))

        summaries.append(
            VesselMatchSummary(
                vessel_name=vessel_name,
                vessel_id=vessel_ids[vessel_name],
                study_a_points_total=len(a),
                study_b_points_total=len(b),
                study_a_points_common_total=total_a,
                study_b_points_common_total=total_b,
                total_common_length_mm=total_len,
                subsections=sub_summaries,
            )
        )
    return common, summaries


# -----------------------------
# VTP writing
# -----------------------------

def fmt_floats(values: Iterable[float]) -> str:
    return " ".join(f"{v:.9g}" for v in values)


def fmt_ints(values: Iterable[int]) -> str:
    return " ".join(str(int(v)) for v in values)


def add_data_array(parent: ET.Element, name: str, values: str, data_type: str = "Float32", ncomp: Optional[int] = None) -> None:
    attrs = {"type": data_type, "Name": name, "format": "ascii"}
    if ncomp is not None:
        attrs["NumberOfComponents"] = str(ncomp)
    arr = ET.SubElement(parent, "DataArray", attrs)
    arr.text = values


def write_polyline_vtp(output_path: Path, line_items: Sequence[Tuple[str, int, int, Sequence[CrossSection]]]) -> None:
    points: List[Tuple[float, float, float]] = []
    point_study_id: List[int] = []
    point_vessel_id: List[int] = []
    point_distance: List[float] = []
    point_original_index: List[int] = []
    connectivity: List[int] = []
    offsets: List[int] = []
    cell_study_id: List[int] = []
    cell_vessel_id: List[int] = []

    for _name, vessel_id, study_id, sections in line_items:
        if len(sections) < 2:
            continue
        start_idx = len(points)
        for s in sections:
            points.append(s.position)
            point_study_id.append(study_id)
            point_vessel_id.append(vessel_id)
            point_distance.append(s.vessel_distance)
            point_original_index.append(s.original_index)
        connectivity.extend(range(start_idx, start_idx + len(sections)))
        offsets.append(len(connectivity))
        cell_study_id.append(study_id)
        cell_vessel_id.append(vessel_id)

    root = ET.Element("VTKFile", {"type": "PolyData", "version": "0.1", "byte_order": "LittleEndian"})
    polydata = ET.SubElement(root, "PolyData")
    piece = ET.SubElement(
        polydata,
        "Piece",
        {
            "NumberOfPoints": str(len(points)),
            "NumberOfVerts": "0",
            "NumberOfLines": str(len(offsets)),
            "NumberOfStrips": "0",
            "NumberOfPolys": "0",
        },
    )
    point_data = ET.SubElement(piece, "PointData", {"Scalars": "vessel_id"})
    add_data_array(point_data, "study_id", fmt_ints(point_study_id), "Int32")
    add_data_array(point_data, "vessel_id", fmt_ints(point_vessel_id), "Int32")
    add_data_array(point_data, "distance_mm", fmt_floats(point_distance), "Float32")
    add_data_array(point_data, "original_index", fmt_ints(point_original_index), "Int32")
    cell_data = ET.SubElement(piece, "CellData", {"Scalars": "vessel_id"})
    add_data_array(cell_data, "study_id", fmt_ints(cell_study_id), "Int32")
    add_data_array(cell_data, "vessel_id", fmt_ints(cell_vessel_id), "Int32")
    points_el = ET.SubElement(piece, "Points")
    add_data_array(points_el, "Points", fmt_floats([c for p in points for c in p]), "Float32", 3)
    lines_el = ET.SubElement(piece, "Lines")
    add_data_array(lines_el, "connectivity", fmt_ints(connectivity), "Int32")
    add_data_array(lines_el, "offsets", fmt_ints(offsets), "Int32")
    tree = ET.ElementTree(root)
    ET.indent(tree, space="  ")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    tree.write(output_path, encoding="utf-8", xml_declaration=True)


def write_common_vtp(output_path: Path, common: Dict[str, List[Tuple[List[CrossSection], List[CrossSection]]]], vessel_ids: Dict[str, int]) -> None:
    line_items: List[Tuple[str, int, int, Sequence[CrossSection]]] = []
    for vessel_name in sorted(common):
        for a_sections, b_sections in common[vessel_name]:
            line_items.append((vessel_name, vessel_ids[vessel_name], 1, a_sections))
            line_items.append((vessel_name, vessel_ids[vessel_name], 2, b_sections))
    write_polyline_vtp(output_path, line_items)


def write_complete_vessels_vtp(output_path: Path, vessels: Dict[str, List[CrossSection]], vessel_ids: Dict[str, int], study_id: int) -> None:
    line_items = [(name, vessel_ids[name], study_id, vessels[name]) for name in sorted(vessels)]
    write_polyline_vtp(output_path, line_items)


# -----------------------------
# SDF processing
# -----------------------------

PLAQUE_COMPONENTS = ("CALC", "LRNC", "NonCALCMATX")
VOLUME_COMPONENTS = ("lumen", "wall", "CALC", "LRNC", "NonCALCMATX")


def require_sitk() -> None:
    if sitk is None:
        raise RuntimeError("SimpleITK is required. Install it with: pip install SimpleITK")


def parse_study_inputs(value: str, label: str) -> StudyInputs:
    parts = [p.strip() for p in str(value).split(",")]
    if len(parts) != 5:
        raise ValueError(
            f"{label} must contain exactly five comma-separated paths: "
            "readings.json,lumen.nrrd,wall.nrrd,plaque.multi.nrrd,wallPartition.multi.nrrd"
        )

    expanded = [expand_path(p) for p in parts]
    return StudyInputs(
        readings=expanded[0],
        lumen=expanded[1],
        wall=expanded[2],
        plaque_multi=expanded[3],
        wall_partition_multi=expanded[4],
    )


def find_multi_nrrd_tool() -> Path:
    here = Path(__file__).resolve().parent
    for name in ("multi_nrrd_reader_writer.py", "multi_nrrd_reader_writer_site_aware.py"):
        p = here / name
        if p.exists():
            return p.resolve()
    raise FileNotFoundError(f"Could not find multi_nrrd_reader_writer.py in {here}")


def run_multi_nrrd_read(input_path: Path, output_dir: Path) -> None:
    input_path = expand_path(input_path)
    output_dir = expand_path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    cmd = [
        sys.executable,
        str(find_multi_nrrd_tool()),
        "read",
        "--input",
        str(input_path),
        "--output-dir",
        str(output_dir),
    ]
    print("Running:", " ".join(cmd))
    subprocess.run(cmd, check=True)


def unpack_study_multi_nrrds(study: StudyInputs, args: argparse.Namespace, study_label: str) -> Dict[str, Path]:
    base = args.inter_out / args.target / study_label
    plaque_dir = base / "plaque"
    wall_partition_dir = base / "wallPartition"
    run_multi_nrrd_read(study.plaque_multi, plaque_dir)
    run_multi_nrrd_read(study.wall_partition_multi, wall_partition_dir)
    return {"base": base, "plaque_dir": plaque_dir, "wall_partition_dir": wall_partition_dir}


def find_matching_nrrd_by_name(directory: Path, name: str) -> Optional[Path]:
    exact = directory / f"{name}.nrrd"
    if exact.exists():
        return exact
    wanted = normalize_name(name).lower()
    for p in sorted(directory.glob("*.nrrd")):
        if normalize_name(p.stem).lower() == wanted:
            return p
    return None


def read_sdf(path: Path) -> Tuple[Any, np.ndarray]:
    require_sitk()
    img = sitk.ReadImage(str(path))
    arr = sitk.GetArrayFromImage(img).astype(np.float32, copy=False)
    return img, arr


def assert_same_geometry_pair(ref_img: Any, ref_arr: np.ndarray, img: Any, arr: np.ndarray, label: str) -> None:
    if arr.shape != ref_arr.shape:
        raise ValueError(f"{label}: shape mismatch {arr.shape} != {ref_arr.shape}")
    if img.GetSize() != ref_img.GetSize():
        raise ValueError(f"{label}: image size mismatch")
    if img.GetSpacing() != ref_img.GetSpacing():
        raise ValueError(f"{label}: spacing mismatch")
    if img.GetOrigin() != ref_img.GetOrigin():
        raise ValueError(f"{label}: origin mismatch")
    if img.GetDirection() != ref_img.GetDirection():
        raise ValueError(f"{label}: direction mismatch")


def voxel_volume_mm3(img: Any) -> float:
    sx, sy, sz = img.GetSpacing()
    return float(sx * sy * sz)


def physical_points_from_zyx_indices(img: Any, zyx_indices: np.ndarray) -> np.ndarray:
    if zyx_indices.size == 0:
        return np.empty((0, 3), dtype=np.float64)
    idx_xyz = zyx_indices[:, [2, 1, 0]].astype(np.float64)
    spacing = np.asarray(img.GetSpacing(), dtype=np.float64)
    origin = np.asarray(img.GetOrigin(), dtype=np.float64)
    direction = np.asarray(img.GetDirection(), dtype=np.float64).reshape(3, 3)
    return origin + (idx_xyz * spacing) @ direction.T


def nearest_centerline_distances(points_xyz: np.ndarray, centerline_xyz: np.ndarray) -> np.ndarray:
    if points_xyz.size == 0:
        return np.empty((0,), dtype=np.float64)
    if cKDTree is not None:
        return cKDTree(centerline_xyz).query(points_xyz, k=1)[0].astype(np.float64, copy=False)
    out = np.empty((points_xyz.shape[0],), dtype=np.float64)
    for start in range(0, points_xyz.shape[0], 25000):
        stop = min(start + 25000, points_xyz.shape[0])
        p = points_xyz[start:stop, None, :]
        c = centerline_xyz[None, :, :]
        out[start:stop] = np.sqrt(np.sum((p - c) ** 2, axis=2)).min(axis=1)
    return out


def normalized(v: np.ndarray) -> np.ndarray:
    n = float(np.linalg.norm(v))
    return v if n < 1e-12 else v / n


def cross_section_normal(section: CrossSection, fallback_tangent: np.ndarray) -> np.ndarray:
    normal = None
    if section.xaxis is not None and section.yaxis is not None:
        x = np.asarray(section.xaxis, dtype=np.float64)
        y = np.asarray(section.yaxis, dtype=np.float64)
        n = np.cross(x, y)
        if np.linalg.norm(n) > 1e-12:
            normal = normalized(n)
    if normal is None:
        normal = normalized(fallback_tangent)
    tangent = normalized(fallback_tangent)
    if np.dot(normal, tangent) < 0:
        normal = -normal
    return normal


def cap_plane_mask(points_xyz: np.ndarray, sections: Sequence[CrossSection], cap_tolerance_mm: float) -> np.ndarray:
    if len(sections) < 2 or points_xyz.size == 0:
        return np.zeros((points_xyz.shape[0],), dtype=bool)
    p_start = np.asarray(sections[0].position, dtype=np.float64)
    p_end = np.asarray(sections[-1].position, dtype=np.float64)
    start_tangent = np.asarray(sections[1].position, dtype=np.float64) - p_start
    end_tangent = p_end - np.asarray(sections[-2].position, dtype=np.float64)
    n_start = cross_section_normal(sections[0], start_tangent)
    n_end = cross_section_normal(sections[-1], end_tangent)
    return ((points_xyz - p_start) @ n_start >= -cap_tolerance_mm) & ((points_xyz - p_end) @ n_end <= cap_tolerance_mm)


def make_common_vessel_mask(
    reference_img: Any,
    vessel_wall_arr: np.ndarray,
    centerline_sections: Sequence[CrossSection],
    sdf_threshold_mm: float,
    centerline_radius_mm: float,
    cap_tolerance_mm: float,
) -> np.ndarray:
    if len(centerline_sections) < 2:
        return np.zeros_like(vessel_wall_arr, dtype=bool)
    candidate_zyx = np.argwhere(vessel_wall_arr <= sdf_threshold_mm)
    if candidate_zyx.size == 0:
        return np.zeros_like(vessel_wall_arr, dtype=bool)
    candidate_xyz = physical_points_from_zyx_indices(reference_img, candidate_zyx)
    cap_keep = cap_plane_mask(candidate_xyz, centerline_sections, cap_tolerance_mm)
    if not np.any(cap_keep):
        return np.zeros_like(vessel_wall_arr, dtype=bool)
    capped_zyx = candidate_zyx[cap_keep]
    capped_xyz = candidate_xyz[cap_keep]
    centerline_xyz = np.asarray([s.position for s in centerline_sections], dtype=np.float64)
    keep = nearest_centerline_distances(capped_xyz, centerline_xyz) <= centerline_radius_mm
    mask = np.zeros_like(vessel_wall_arr, dtype=bool)
    if np.any(keep):
        kept = capped_zyx[keep]
        mask[kept[:, 0], kept[:, 1], kept[:, 2]] = True
    return mask


def safe_name(name: str) -> str:
    return re.sub(r"[^A-Za-z0-9_]+", "_", str(name)).strip("_") or "Unknown"



def combined_label_value(vessel_index: int, subsection_index: int, component_index: int) -> int:
    """
    Encode labels into one combined per-study QA label map.

    Encoding:
      label = vessel_index * 1000 + subsection_index * 10 + component_index

    component_index:
      1 lumen
      2 wall
      3 CALC
      4 LRNC
      5 NonCALCMATX
    """
    return int(vessel_index * 1000 + subsection_index * 10 + component_index)


COMPONENT_LABEL_IDS = {
    "lumen": 1,
    "wall": 2,
    "CALC": 3,
    "LRNC": 4,
    "NonCALCMATX": 5,
}


def write_combined_label_map(combined_label: np.ndarray, reference_img: Any, output_path: Path) -> Dict[str, str]:
    require_sitk()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    label_img = sitk.GetImageFromArray(combined_label.astype(np.uint16))
    label_img.CopyInformation(reference_img)
    sitk.WriteImage(label_img, str(output_path))
    return {"combined_label_map_nii_gz": str(output_path)}


def write_mask_outputs(
    mask: np.ndarray,
    reference_img: Any,
    output_stem: Path,
    write_label: bool = False,
) -> Dict[str, str]:
    """
    Write restricted component SDF NRRD. Individual binary label maps are
    optional; by default the pipeline writes one combined per-study label map.
    """
    require_sitk()

    label_img = sitk.GetImageFromArray(mask.astype(np.uint8))
    label_img.CopyInformation(reference_img)
    sdf_img = sitk.SignedMaurerDistanceMap(
        label_img,
        insideIsPositive=False,
        squaredDistance=False,
        useImageSpacing=True,
    )

    sdf_path = output_stem.with_suffix(".sdf.nrrd")
    sdf_path.parent.mkdir(parents=True, exist_ok=True)
    sitk.WriteImage(sdf_img, str(sdf_path))

    out = {"sdf_nrrd": str(sdf_path)}
    if write_label:
        label_path = output_stem.with_suffix(".nii.gz")
        sitk.WriteImage(label_img, str(label_path))
        out["label_map_nii_gz"] = str(label_path)

    return out



def require_vtk() -> None:
    if vtk is None or vtk_numpy_support is None:
        raise RuntimeError(
            "VTK is required for --volume-mode subvoxel. "
            "Install it with: pip install vtk"
        )


def image_geometry_to_vtk_image(reference_img: Any, arr_zyx: np.ndarray) -> Any:
    require_vtk()

    arr_xyz = np.transpose(arr_zyx.astype(np.float32, copy=False), (2, 1, 0))
    arr_flat = np.ascontiguousarray(arr_xyz).ravel(order="F")

    vtk_arr = vtk_numpy_support.numpy_to_vtk(
        num_array=arr_flat,
        deep=True,
        array_type=vtk.VTK_FLOAT,
    )
    vtk_arr.SetName("sdf")

    image = vtk.vtkImageData()
    image.SetDimensions(reference_img.GetSize())
    image.SetSpacing(reference_img.GetSpacing())
    image.SetOrigin(reference_img.GetOrigin())
    image.GetPointData().SetScalars(vtk_arr)
    return image


def subvoxel_volume_from_sdf_arr(sdf_arr: np.ndarray, reference_img: Any) -> float:
    """
    Compute sub-voxel volume from the SDF zero level set using:
      vtkFlyingEdges3D -> vtkTriangleFilter -> vtkCleanPolyData -> vtkMassProperties

    SDF convention:
      negative inside, positive outside.
    """
    require_vtk()

    if not np.any(sdf_arr <= 0.0):
        return 0.0

    vtk_image = image_geometry_to_vtk_image(reference_img, sdf_arr)

    surface = vtk.vtkFlyingEdges3D()
    surface.SetInputData(vtk_image)
    surface.SetValue(0, 0.0)
    surface.ComputeNormalsOn()
    surface.Update()

    poly = surface.GetOutput()
    if poly is None or poly.GetNumberOfPoints() == 0 or poly.GetNumberOfCells() == 0:
        return 0.0

    tri = vtk.vtkTriangleFilter()
    tri.SetInputData(poly)
    tri.Update()

    clean = vtk.vtkCleanPolyData()
    clean.SetInputData(tri.GetOutput())
    clean.Update()

    mass = vtk.vtkMassProperties()
    mass.SetInputData(clean.GetOutput())
    mass.Update()

    volume = float(mass.GetVolume())
    if not np.isfinite(volume):
        return 0.0
    return abs(volume)


def signed_distance_to_cap_planes_arr(
    reference_img: Any,
    shape_zyx: Tuple[int, int, int],
    sections: Sequence[CrossSection],
    cap_tolerance_mm: float,
) -> np.ndarray:
    if len(sections) < 2:
        return np.ones(shape_zyx, dtype=np.float32)

    all_zyx = np.indices(shape_zyx, dtype=np.int32).reshape(3, -1).T
    points_xyz = physical_points_from_zyx_indices(reference_img, all_zyx)

    start = sections[0]
    end = sections[-1]
    p_start = np.asarray(start.position, dtype=np.float64)
    p_end = np.asarray(end.position, dtype=np.float64)

    start_tangent = np.asarray(sections[1].position, dtype=np.float64) - p_start
    end_tangent = p_end - np.asarray(sections[-2].position, dtype=np.float64)

    n_start = cross_section_normal(start, start_tangent)
    n_end = cross_section_normal(end, end_tangent)

    proximal_phi = -((points_xyz - p_start) @ n_start) - cap_tolerance_mm
    distal_phi = ((points_xyz - p_end) @ n_end) - cap_tolerance_mm

    return np.maximum(proximal_phi, distal_phi).reshape(shape_zyx).astype(np.float32)


def centerline_radius_sdf_arr(
    reference_img: Any,
    shape_zyx: Tuple[int, int, int],
    sections: Sequence[CrossSection],
    centerline_radius_mm: float,
) -> np.ndarray:
    if len(sections) < 2:
        return np.ones(shape_zyx, dtype=np.float32)

    all_zyx = np.indices(shape_zyx, dtype=np.int32).reshape(3, -1).T
    points_xyz = physical_points_from_zyx_indices(reference_img, all_zyx)
    centerline_xyz = np.asarray([s.position for s in sections], dtype=np.float64)

    nearest_dist = nearest_centerline_distances(points_xyz, centerline_xyz)
    return (nearest_dist - float(centerline_radius_mm)).reshape(shape_zyx).astype(np.float32)


def make_common_section_sdf(
    reference_img: Any,
    vessel_wall_arr: np.ndarray,
    sections: Sequence[CrossSection],
    centerline_radius_mm: float,
    cap_tolerance_mm: float,
) -> np.ndarray:
    """
    Level-set intersection for the closed common vessel section:
      max(vessel_wall_sdf, cap_planes_sdf, centerline_radius_sdf)
    """
    cap_phi = signed_distance_to_cap_planes_arr(
        reference_img,
        vessel_wall_arr.shape,
        sections,
        cap_tolerance_mm,
    )
    radius_phi = centerline_radius_sdf_arr(
        reference_img,
        vessel_wall_arr.shape,
        sections,
        centerline_radius_mm,
    )
    return np.maximum.reduce([
        vessel_wall_arr.astype(np.float32, copy=False),
        cap_phi,
        radius_phi,
    ])


def component_common_sdf_arr(component_arr: np.ndarray, section_sdf: np.ndarray) -> np.ndarray:
    """
    Level-set intersection:
      component inside common vessel section = max(component_sdf, section_sdf)
    """
    return np.maximum(
        component_arr.astype(np.float32, copy=False),
        section_sdf.astype(np.float32, copy=False),
    )




def write_sdf_array_output(
    sdf_arr: np.ndarray,
    reference_img: Any,
    output_stem: Path,
) -> Dict[str, str]:
    """
    Write an SDF array directly without converting through a binary label map.

    This preserves sub-voxel smoothness of original SDF level-set operations.
    Use this for --volume-mode subvoxel QA SDF outputs.
    """
    require_sitk()
    sdf_path = output_stem.with_suffix(".sdf.nrrd")
    sdf_path.parent.mkdir(parents=True, exist_ok=True)

    sdf_img = sitk.GetImageFromArray(sdf_arr.astype(np.float32, copy=False))
    sdf_img.CopyInformation(reference_img)
    sitk.WriteImage(sdf_img, str(sdf_path))

    return {"sdf_nrrd": str(sdf_path)}


def compute_component_volume(
    component_arr: np.ndarray,
    section_mask: np.ndarray,
    threshold: float,
    voxel_volume: float,
    volume_mode: str = "voxel",
    component_common_sdf: Optional[np.ndarray] = None,
    reference_img: Optional[Any] = None,
) -> Tuple[np.ndarray, int, float]:
    mask = (component_arr <= threshold) & section_mask
    count = int(np.count_nonzero(mask))

    if volume_mode == "voxel":
        return mask, count, float(count * voxel_volume)

    if volume_mode == "subvoxel":
        if component_common_sdf is None or reference_img is None:
            raise ValueError("component_common_sdf and reference_img are required for subvoxel volume.")
        return mask, count, subvoxel_volume_from_sdf_arr(component_common_sdf, reference_img)

    raise ValueError(f"Unsupported volume mode: {volume_mode}")


def process_study_sdf_volumes(
    args: argparse.Namespace,
    study: StudyInputs,
    study_label: str,
    study_id: int,
    unpacked: Dict[str, Path],
    common: Dict[str, List[Tuple[List[CrossSection], List[CrossSection]]]],
) -> Tuple[Dict[str, Any], Dict[str, str]]:
    if args.volume_mode == "subvoxel":
        require_vtk()

    lumen_img, lumen_arr = read_sdf(study.lumen)
    wall_img, wall_arr = read_sdf(study.wall)
    assert_same_geometry_pair(lumen_img, lumen_arr, wall_img, wall_arr, f"{study_label} wall")

    plaque_arrays: Dict[str, Tuple[Any, np.ndarray]] = {}
    for component in PLAQUE_COMPONENTS:
        p = unpacked["plaque_dir"] / f"{component}.nrrd"
        if p.exists():
            img, arr = read_sdf(p)
            assert_same_geometry_pair(lumen_img, lumen_arr, img, arr, f"{study_label} {component}")
            plaque_arrays[component] = (img, arr)
        else:
            print(f"Warning: missing plaque component {component}: {p}")

    voxel_vol = voxel_volume_mm3(lumen_img)
    volume_summary: Dict[str, Any] = {}
    generated: Dict[str, str] = {}

    combined_label = np.zeros_like(lumen_arr, dtype=np.uint16)
    combined_label_legend: Dict[str, Any] = {
        "encoding": "label = vessel_index * 1000 + subsection_index * 10 + component_index",
        "component_label_ids": COMPONENT_LABEL_IDS,
        "labels": {},
    }

    for vessel_index, vessel_name in enumerate(sorted(common), start=1):
        vessel_wall_path = find_matching_nrrd_by_name(unpacked["wall_partition_dir"], vessel_name)
        if vessel_wall_path is None:
            print(f"Warning: no wallPartition vessel SDF found for {vessel_name}; skipping")
            continue

        vessel_wall_img, vessel_wall_arr = read_sdf(vessel_wall_path)
        assert_same_geometry_pair(
            lumen_img,
            lumen_arr,
            vessel_wall_img,
            vessel_wall_arr,
            f"{study_label} wallPartition {vessel_name}",
        )

        vessel_key = safe_name(vessel_name)
        volume_summary[vessel_name] = {
            "wall_partition_sdf": str(vessel_wall_path),
            "subsections": [],
            "total": {
                "common_length_mm": 0.0,
                "components": {},
            },
        }

        total_masks: Dict[str, np.ndarray] = {
            component: np.zeros_like(lumen_arr, dtype=bool)
            for component in VOLUME_COMPONENTS
        }
        # For subvoxel totals, preserve level-set union:
        # union(phi_1, phi_2, ...) = min(phi_1, phi_2, ...)
        total_sdfs: Dict[str, Optional[np.ndarray]] = {
            component: None for component in VOLUME_COMPONENTS
        }

        for subsection_idx, (sections_a, sections_b) in enumerate(common[vessel_name], start=1):
            sections = sections_a if study_id == 1 else sections_b
            if len(sections) < 2:
                continue

            subsection_length = max(0.0, sections[-1].vessel_distance - sections[0].vessel_distance)
            section_mask = make_common_vessel_mask(
                reference_img=lumen_img,
                vessel_wall_arr=vessel_wall_arr,
                centerline_sections=sections,
                sdf_threshold_mm=args.sdf_threshold_mm,
                centerline_radius_mm=args.centerline_radius_mm,
                cap_tolerance_mm=args.cap_tolerance_mm,
            )

            section_sdf = make_common_section_sdf(
                reference_img=lumen_img,
                vessel_wall_arr=vessel_wall_arr,
                sections=sections,
                centerline_radius_mm=args.centerline_radius_mm,
                cap_tolerance_mm=args.cap_tolerance_mm,
            )

            subsection_dir = (
                args.inter_out
                / args.target
                / study_label
                / "common_sections"
                / vessel_key
                / f"section_{subsection_idx:03d}"
            )

            payload = {
                "subsection_id": subsection_idx,
                "common_start_mm": float(sections[0].vessel_distance),
                "common_end_mm": float(sections[-1].vessel_distance),
                "common_length_mm": float(subsection_length),
                "common_vessel_mask_voxel_count": int(np.count_nonzero(section_mask)),
                "common_vessel_mask_volume_mm3": float(np.count_nonzero(section_mask) * voxel_vol),
                "components": {},
            }

            sources: Dict[str, Tuple[Any, np.ndarray]] = {
                "lumen": (lumen_img, lumen_arr),
                "wall": (wall_img, wall_arr),
            }
            sources.update(plaque_arrays)

            for component, (component_img, component_arr) in sources.items():
                component_common_sdf = component_common_sdf_arr(component_arr, section_sdf)

                component_mask, count, vol = compute_component_volume(
                    component_arr,
                    section_mask,
                    args.sdf_threshold_mm,
                    voxel_vol,
                    volume_mode=args.volume_mode,
                    component_common_sdf=component_common_sdf,
                    reference_img=component_img,
                )
                total_masks[component] |= component_mask

                if args.volume_mode == "subvoxel":
                    if total_sdfs[component] is None:
                        total_sdfs[component] = component_common_sdf.copy()
                    else:
                        total_sdfs[component] = np.minimum(total_sdfs[component], component_common_sdf)

                if args.volume_mode == "subvoxel":
                    # Preserve the level-set SDF generated by intersection:
                    # max(component_sdf, common_section_sdf).
                    output_paths = write_sdf_array_output(
                        component_common_sdf,
                        component_img,
                        subsection_dir / f"{safe_name(component)}_common",
                    )
                else:
                    # Voxel mode preserves existing behavior: binary mask -> SDF.
                    output_paths = write_mask_outputs(
                        component_mask,
                        component_img,
                        subsection_dir / f"{safe_name(component)}_common",
                        write_label=False,
                    )

                component_label_id = COMPONENT_LABEL_IDS[component]
                label_value = combined_label_value(
                    vessel_index=vessel_index,
                    subsection_index=subsection_idx,
                    component_index=component_label_id,
                )
                combined_label[component_mask] = label_value
                combined_label_legend["labels"][str(label_value)] = {
                    "study": study_label,
                    "vessel": vessel_name,
                    "vessel_index": vessel_index,
                    "subsection_id": subsection_idx,
                    "component": component,
                    "component_index": component_label_id,
                }

                payload["components"][component] = {
                    "voxel_count": count,
                    "volume_mm3": vol,
                    "volume_method": args.volume_mode,
                    "combined_label_value": label_value,
                    "outputs": output_paths,
                }

                generated[
                    f"{study_label}_{vessel_key}_section_{subsection_idx:03d}_{component}_sdf"
                ] = output_paths["sdf_nrrd"]

            volume_summary[vessel_name]["subsections"].append(payload)
            volume_summary[vessel_name]["total"]["common_length_mm"] += float(subsection_length)

        total_dir = args.inter_out / args.target / study_label / "common_sections" / vessel_key / "total"
        for component, mask in total_masks.items():
            if component not in ("lumen", "wall") and component not in plaque_arrays:
                continue

            ref_img = lumen_img if component == "lumen" else wall_img
            if component in plaque_arrays:
                ref_img = plaque_arrays[component][0]

            count = int(np.count_nonzero(mask))
            vol = float(count * voxel_vol)
            if args.volume_mode == "subvoxel" and total_sdfs.get(component) is not None:
                total_sdf_arr = total_sdfs[component]
                output_paths = write_sdf_array_output(
                    total_sdf_arr,
                    ref_img,
                    total_dir / f"{safe_name(component)}_common_total",
                )
                total_volume = subvoxel_volume_from_sdf_arr(total_sdf_arr, ref_img)
            else:
                output_paths = write_mask_outputs(
                    mask,
                    ref_img,
                    total_dir / f"{safe_name(component)}_common_total",
                    write_label=False,
                )
                total_volume = vol

            volume_summary[vessel_name]["total"]["components"][component] = {
                "voxel_count": count,
                "volume_mm3": total_volume,
                "volume_method": args.volume_mode,
                "outputs": output_paths,
            }

            generated[f"{study_label}_{vessel_key}_total_{component}_sdf"] = output_paths["sdf_nrrd"]

    combined_label_path = (
        args.inter_out
        / args.target
        / study_label
        / "common_sections"
        / f"{study_label}_combined_common_sections_label.nii.gz"
    )
    combined_label_outputs = write_combined_label_map(combined_label, lumen_img, combined_label_path)
    generated[f"{study_label}_combined_label_map"] = combined_label_outputs["combined_label_map_nii_gz"]

    volume_summary["_combined_label_map"] = {
        "outputs": combined_label_outputs,
        "legend": combined_label_legend,
    }

    return volume_summary, generated


def process_sdf_volume_pipeline(
    args: argparse.Namespace,
    study_a: StudyInputs,
    study_b: StudyInputs,
    common: Dict[str, List[Tuple[List[CrossSection], List[CrossSection]]]],
) -> Tuple[Dict[str, Any], Dict[str, str]]:
    unpacked_a = unpack_study_multi_nrrds(study_a, args, "study_a")
    unpacked_b = unpack_study_multi_nrrds(study_b, args, "study_b")
    volume_results: Dict[str, Any] = {}
    generated: Dict[str, str] = {
        "study_a_plaque_unpacked_dir": str(unpacked_a["plaque_dir"]),
        "study_a_wall_partition_unpacked_dir": str(unpacked_a["wall_partition_dir"]),
        "study_b_plaque_unpacked_dir": str(unpacked_b["plaque_dir"]),
        "study_b_wall_partition_unpacked_dir": str(unpacked_b["wall_partition_dir"]),
    }
    volume_results["study_a"], files_a = process_study_sdf_volumes(args, study_a, "study_a", 1, unpacked_a, common)
    volume_results["study_b"], files_b = process_study_sdf_volumes(args, study_b, "study_b", 2, unpacked_b, common)
    generated.update(files_a)
    generated.update(files_b)
    return volume_results, generated


# -----------------------------
# Summary / CLI
# -----------------------------


def get_component_volume_from_study(
    volume_results: Dict[str, Any],
    study_key: str,
    vessel_name: str,
    component: str,
    subsection_index: Optional[int] = None,
) -> float:
    study_payload = volume_results.get(study_key, {})
    vessel_payload = study_payload.get(vessel_name, {})

    if subsection_index is None:
        return float(
            vessel_payload
            .get("total", {})
            .get("components", {})
            .get(component, {})
            .get("volume_mm3", 0.0)
        )

    subsections = vessel_payload.get("subsections", [])
    if subsection_index < 1 or subsection_index > len(subsections):
        return 0.0

    return float(
        subsections[subsection_index - 1]
        .get("components", {})
        .get(component, {})
        .get("volume_mm3", 0.0)
    )


def percentage_delta(numerator_value: float, denominator_value: float) -> Optional[float]:
    """
    Percentage delta of numerator relative to denominator:
      ((numerator - denominator) / denominator) * 100

    Returns None when denominator is effectively zero.
    """
    if abs(denominator_value) < 1e-12:
        return None
    return ((numerator_value - denominator_value) / denominator_value) * 100.0


def format_percent(value: Optional[float]) -> str:
    if value is None:
        return "NA"
    return f"{value:.2f}%"


def format_volume_row(name: str, study_a_value: float, study_b_value: float) -> str:
    diff = study_a_value - study_b_value
    return f"{name:<32}{study_a_value:>12.2f}{study_b_value:>12.2f}{diff:>12.2f}"


def format_total_volume_row(name: str, study_a_value: float, study_b_value: float) -> str:
    diff = study_a_value - study_b_value
    pct_a_vs_b = percentage_delta(study_a_value, study_b_value)
    pct_b_vs_a = percentage_delta(study_b_value, study_a_value)

    return (
        f"{name:<32}"
        f"{study_a_value:>12.2f}"
        f"{study_b_value:>12.2f}"
        f"{diff:>12.2f}"
        f"{format_percent(pct_a_vs_b):>14}"
        f"{format_percent(pct_b_vs_a):>14}"
    )


def write_volume_comparison_txt(
    output_path: Path,
    volume_results: Dict[str, Any],
    summaries: Sequence[VesselMatchSummary],
) -> None:
    """
    Write side-by-side volume comparison table.

    Units are mm3 because volumes are computed as:
      voxel_count * spacing_x * spacing_y * spacing_z

    For component total rows only, include:
      % Delta A vs B = ((Study A - Study B) / Study B) * 100
      % Delta B vs A = ((Study B - Study A) / Study A) * 100
    """
    output_path.parent.mkdir(parents=True, exist_ok=True)

    sections = [
        ("Wall", "wall"),
        ("Lumen", "lumen"),
        ("Plaque - CALC", "CALC"),
        ("Plaque - LRNC", "LRNC"),
        ("Plaque - NonCALCMATX", "NonCALCMATX"),
    ]

    vessel_summaries = sorted(summaries, key=lambda s: s.vessel_name)
    line_width = 96
    lines: List[str] = []

    lines.append("Volume comparison for common vessel sections")
    lines.append("Units: mm^3")
    lines.append("Percent deltas are shown only for component total rows.")
    lines.append("% Delta A vs B = ((Study A - Study B) / Study B) * 100")
    lines.append("% Delta B vs A = ((Study B - Study A) / Study A) * 100")
    lines.append("")

    for title, component in sections:
        lines.append(
            f"{'Vessel':<32}"
            f"{'Study A':>12}"
            f"{'Study B':>12}"
            f"{'Diff (A-B)':>12}"
            f"{'%Delta A/B':>14}"
            f"{'%Delta B/A':>14}"
        )
        lines.append("-" * line_width)
        lines.append(title)

        total_a = 0.0
        total_b = 0.0

        for vessel_summary in vessel_summaries:
            vessel_name = vessel_summary.vessel_name

            vessel_a = get_component_volume_from_study(
                volume_results,
                "study_a",
                vessel_name,
                component,
                subsection_index=None,
            )
            vessel_b = get_component_volume_from_study(
                volume_results,
                "study_b",
                vessel_name,
                component,
                subsection_index=None,
            )

            total_a += vessel_a
            total_b += vessel_b

            lines.append(format_volume_row(vessel_name, vessel_a, vessel_b))

            for subsection in vessel_summary.subsections:
                section_a = get_component_volume_from_study(
                    volume_results,
                    "study_a",
                    vessel_name,
                    component,
                    subsection_index=subsection.subsection_id,
                )
                section_b = get_component_volume_from_study(
                    volume_results,
                    "study_b",
                    vessel_name,
                    component,
                    subsection_index=subsection.subsection_id,
                )
                lines.append(
                    format_volume_row(
                        f"    section {subsection.subsection_id}",
                        section_a,
                        section_b,
                    )
                )

            lines.append("")

        if title.startswith("Plaque - "):
            total_label = f"Total - {component}"
        else:
            total_label = f"Total - {component}"

        lines.append(format_total_volume_row(total_label, total_a, total_b))
        lines.append("-" * line_width)
        lines.append("")

    output_path.write_text("\n".join(lines), encoding="utf-8")


def format_xyz(point: Tuple[float, float, float]) -> str:
    return f"({point[0]:.3f}, {point[1]:.3f}, {point[2]:.3f})"


def find_common_pair_sections(
    common: Dict[str, List[Tuple[List[CrossSection], List[CrossSection]]]],
    vessel_name: str,
    subsection_id: int,
) -> Optional[Tuple[List[CrossSection], List[CrossSection]]]:
    subsections = common.get(vessel_name, [])
    if subsection_id < 1 or subsection_id > len(subsections):
        return None
    return subsections[subsection_id - 1]


def write_common_sections_txt(
    output_path: Path,
    summaries: Sequence[VesselMatchSummary],
    common: Dict[str, List[Tuple[List[CrossSection], List[CrossSection]]]],
) -> None:
    """
    Write common cross-section geometry report.

    Includes:
      - vessel aggregate total common length
      - per-subsection common length
      - Study A start/end XYZ coordinates
      - Study B start/end XYZ coordinates
    """
    output_path.parent.mkdir(parents=True, exist_ok=True)

    line_width = 156
    lines: List[str] = []

    lines.append("Common vessel cross-section geometry summary")
    lines.append("Units: coordinates in image physical space; lengths in mm")
    lines.append("")
    lines.append(
        f"{'Vessel / Section':<34}"
        f"{'Length':>10}"
        f"{'A Start XYZ':>28}"
        f"{'A End XYZ':>28}"
        f"{'B Start XYZ':>28}"
        f"{'B End XYZ':>28}"
    )
    lines.append("-" * line_width)

    grand_total = 0.0

    for vessel_summary in sorted(summaries, key=lambda s: s.vessel_name):
        vessel_name = vessel_summary.vessel_name
        grand_total += float(vessel_summary.total_common_length_mm)

        # Aggregate vessel row: show aggregate length and leave coordinates blank
        # because multiple discrete subsections may have multiple start/end caps.
        lines.append(
            f"{vessel_name:<34}"
            f"{vessel_summary.total_common_length_mm:>10.2f}"
            f"{'':>28}"
            f"{'':>28}"
            f"{'':>28}"
            f"{'':>28}"
        )

        for subsection in vessel_summary.subsections:
            pair = find_common_pair_sections(common, vessel_name, subsection.subsection_id)
            if pair is None:
                continue

            a_sections, b_sections = pair
            if not a_sections or not b_sections:
                continue

            a_start = a_sections[0].position
            a_end = a_sections[-1].position
            b_start = b_sections[0].position
            b_end = b_sections[-1].position

            lines.append(
                f"{'    section ' + str(subsection.subsection_id):<34}"
                f"{subsection.common_length_mm:>10.2f}"
                f"{format_xyz(a_start):>28}"
                f"{format_xyz(a_end):>28}"
                f"{format_xyz(b_start):>28}"
                f"{format_xyz(b_end):>28}"
            )

        lines.append("")

    lines.append("-" * line_width)
    lines.append(
        f"{'Total common length':<34}"
        f"{grand_total:>10.2f}"
        f"{'':>28}"
        f"{'':>28}"
        f"{'':>28}"
        f"{'':>28}"
    )
    lines.append("-" * line_width)

    output_path.write_text("\n".join(lines), encoding="utf-8")


def write_summary(
    output_path: Path,
    summaries: Sequence[VesselMatchSummary],
    args: argparse.Namespace,
    generated_files: Dict[str, str],
    vessel_ids: Dict[str, int],
    study_a: StudyInputs,
    study_b: StudyInputs,
    volume_results: Optional[Dict[str, Any]] = None,
) -> None:
    payload = {
        "target": args.target,
        "volume_mode": args.volume_mode,
        "study_a": {k: str(v) for k, v in asdict(study_a).items()},
        "study_b": {k: str(v) for k, v in asdict(study_b).items()},
        "out": str(args.out),
        "inter_out": str(args.inter_out),
        "tolerance_mm": args.tolerance_mm,
        "min_common_length_mm": args.min_common_length_mm,
        "max_common_gap_mm": args.max_common_gap_mm,
        "centerline_radius_mm": args.centerline_radius_mm,
        "cap_tolerance_mm": args.cap_tolerance_mm,
        "sdf_threshold_mm": args.sdf_threshold_mm,
        "vessel_ids": vessel_ids,
        "generated_files": generated_files,
        "matched_vessels": [asdict(s) for s in summaries],
        "volume_results": volume_results or {},
    }
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2)


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Find common vessel sections and compute SDF volumes for two studies.")
    parser.add_argument("--study-a", required=True, help="readings.json,lumen.nrrd,wall.nrrd,plaque.multi.nrrd,wallPartition.multi.nrrd")
    parser.add_argument("--study-b", required=True, help="readings.json,lumen.nrrd,wall.nrrd,plaque.multi.nrrd,wallPartition.multi.nrrd")
    parser.add_argument("--target", required=True, choices=["left", "right"], help="Target side/tree: left or right")
    parser.add_argument("--out", required=True, type=Path, help="Final output directory")
    parser.add_argument("--inter_out", required=True, type=Path, help="Intermediate/debug output directory")
    parser.add_argument(
        "--volume-mode",
        choices=["voxel", "subvoxel"],
        default="voxel",
        help=(
            "Volume calculation mode. 'voxel' preserves current behavior using "
            "voxel_count * voxel_volume. 'subvoxel' extracts the SDF zero-level "
            "surface using vtkFlyingEdges3D and computes mesh volume with vtkMassProperties."
        ),
    )
    parser.add_argument("--sdf-threshold-mm", type=float, default=0.0, help="SDF inside threshold; default sdf <= 0")
    parser.add_argument("--centerline-radius-mm", type=float, default=25.0, help="Safety radius around centerline samples")
    parser.add_argument("--cap-tolerance-mm", type=float, default=0.5, help="Tolerance for proximal/distal cap planes")
    parser.add_argument("--tolerance-mm", type=float, default=0.75, help="Centerline common-distance match tolerance")
    parser.add_argument("--min-common-length-mm", type=float, default=5.0, help="Minimum accepted discrete common subsection length")
    parser.add_argument("--max-common-gap-mm", type=float, default=None, help="Maximum gap within a common subsection; default derives from sampling")
    return parser


def main() -> int:
    parser = build_arg_parser()
    args = parser.parse_args()
    if args.tolerance_mm < 0:
        parser.error("--tolerance-mm must be non-negative")
    if args.min_common_length_mm < 0:
        parser.error("--min-common-length-mm must be non-negative")
    if args.max_common_gap_mm is not None and args.max_common_gap_mm < 0:
        parser.error("--max-common-gap-mm must be non-negative")

    args.out = expand_path(args.out)
    args.inter_out = expand_path(args.inter_out)

    args.out = expand_path(args.out)
    args.inter_out = expand_path(args.inter_out)

    study_a = parse_study_inputs(args.study_a, "--study-a")
    study_b = parse_study_inputs(args.study_b, "--study-b")

    vessels_a = parse_readings_json(study_a.readings)
    vessels_b = parse_readings_json(study_b.readings)
    vessel_ids = build_vessel_ids(vessels_a, vessels_b)
    common, summaries = compute_common_sections(
        vessels_a,
        vessels_b,
        args.tolerance_mm,
        vessel_ids,
        args.min_common_length_mm,
        args.max_common_gap_mm,
    )

    args.out.mkdir(parents=True, exist_ok=True)
    target_inter = args.inter_out / args.target
    target_inter.mkdir(parents=True, exist_ok=True)
    common_vtp = target_inter / "common_vessel_sections.vtp"
    study_a_vtp = target_inter / "study_a_complete_vessels.vtp"
    study_b_vtp = target_inter / "study_b_complete_vessels.vtp"
    summary_json = args.out / "common_vessel_sections_summary.json"
    volume_txt = args.out / "common_vessel_sections_volume_summary.txt"
    common_sections_txt = args.out / "common_vessel_sections_geometry_summary.txt"

    write_common_vtp(common_vtp, common, vessel_ids)
    write_complete_vessels_vtp(study_a_vtp, vessels_a, vessel_ids, study_id=1)
    write_complete_vessels_vtp(study_b_vtp, vessels_b, vessel_ids, study_id=2)

    volume_results, sdf_generated_files = process_sdf_volume_pipeline(args, study_a, study_b, common)
    generated_files = {
        "common_vessel_sections_vtp": str(common_vtp),
        "study_a_complete_vessels_vtp": str(study_a_vtp),
        "study_b_complete_vessels_vtp": str(study_b_vtp),
        "summary_json": str(summary_json),
        "volume_summary_txt": str(volume_txt),
        "common_sections_geometry_txt": str(common_sections_txt),
    }
    generated_files.update(sdf_generated_files)

    write_volume_comparison_txt(volume_txt, volume_results, summaries)
    write_common_sections_txt(common_sections_txt, summaries, common)

    write_summary(summary_json, summaries, args, generated_files, vessel_ids, study_a, study_b, volume_results)

    print(f"Wrote common sections VTP: {common_vtp}")
    print(f"Wrote complete study A VTP: {study_a_vtp}")
    print(f"Wrote complete study B VTP: {study_b_vtp}")
    print(f"Wrote summary: {summary_json}")
    print(f"Wrote volume TXT summary: {volume_txt}")
    print(f"Wrote common sections geometry TXT summary: {common_sections_txt}")
    print(f"Volume mode: {args.volume_mode}")
    print("Computed SDF volumes and wrote restricted component label/SDF outputs.")
    print(f"Matched vessels: {len(summaries)}")
    for s in summaries:
        print(f"  vessel_id={s.vessel_id} {s.vessel_name}: {s.total_common_length_mm:.2f} mm across {len(s.subsections)} subsection(s)")
        for sub in s.subsections:
            print(f"    section {sub.subsection_id}: {sub.common_length_mm:.2f} mm, A points={sub.study_a_points_common}, B points={sub.study_b_points_common}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
