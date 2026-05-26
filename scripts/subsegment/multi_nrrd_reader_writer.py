"""
Multi NRRD reader/writer utility for composition.multi.nrrd files.

Fix:
  Extracted 3D NRRDs are written with valid 3D metadata:
    dimension: 3
    space dimension: 3
    space directions: 3 vectors, each with 3 components
    space origin: 3 components

This avoids Slicer errors caused by preserving 4D vectors such as:
  (0.390625,0,0,0)
inside a 3D NRRD file.

Expected region order when slice site metadata is absent:
  CALC, LRNC, IPH, PVAT, MATX, FIBL, NonCALCMATX

Site-aware behavior:
  If a multi-NRRD contains slice0_site, slice1_site, ... metadata, read mode
  writes individual files named by site instead of region.

  In write mode, site metadata is written when provided through --sites or
  --slice-metadata-json. When site names are used, input files are read by
  site name and output order follows the provided list/JSON order.

Requires:
  pip install pynrrd
"""

import argparse
import json
import re
import sys
from pathlib import Path

import numpy as np
import nrrd


REGIONS = [
    "CALC",
    "LRNC",
    "IPH",
    "PVAT",
    "MATX",
    "FIBL",
    "NonCALCMATX",
]


def sanitize_filename(name):
    safe = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(name).strip())
    if not safe:
        raise ValueError(f"Invalid empty region/site name: {name!r}")
    return safe


def get_custom_header_value(header, key, default=None):
    for candidate in (key, key + ":", key + ":="):
        if candidate in header:
            return header[candidate]
    return default


def parse_slice_metadata(header, number_of_slices):
    """
    Parse per-slice metadata from a multi-NRRD header.

    Backward compatible behavior:
      - If slice{i}_site is absent, use slice{i}_region / default REGIONS.

    New site-aware behavior:
      - If any slice{i}_site is present, that is treated as the primary slice
        name for file naming.
      - Region metadata is still parsed and preserved when available.
    """
    any_site_present = any(
        get_custom_header_value(header, f"slice{i}_site", None) is not None
        for i in range(number_of_slices)
    )

    metadata = []
    for i in range(number_of_slices):
        region = get_custom_header_value(header, f"slice{i}_region", None)
        site = get_custom_header_value(header, f"slice{i}_site", None)
        pixel = get_custom_header_value(
            header,
            f"slice{i}_pixelInterpretation",
            "SignedDistanceFunction",
        )
        target = get_custom_header_value(header, f"slice{i}_targetPath", "")

        if region is None:
            region = REGIONS[i] if i < len(REGIONS) else f"slice{i}"

        primary_name = site if any_site_present and site is not None else region

        metadata.append(
            {
                "region": str(region),
                "site": None if site is None else str(site),
                "primaryName": str(primary_name),
                "usesSiteName": bool(any_site_present and site is not None),
                "pixelInterpretation": str(pixel),
                "targetPath": str(target),
            }
        )
    return metadata


def as_float_list(value):
    return [float(v) for v in list(value)]


def get_3d_header_from_4d_header(header4d):
    header3d = {}

    header3d["type"] = header4d.get("type", "float")
    header3d["dimension"] = 3
    header3d["space dimension"] = 3

    if "space directions" not in header4d:
        raise ValueError("Input multi-NRRD is missing 'space directions'")

    dirs4 = header4d["space directions"]
    dirs3 = []
    for d in range(3):
        v = as_float_list(dirs4[d])
        if len(v) < 3:
            raise ValueError(f"Invalid 4D direction vector: {v}")
        dirs3.append([v[0], v[1], v[2]])
    header3d["space directions"] = dirs3

    if "space origin" not in header4d:
        raise ValueError("Input multi-NRRD is missing 'space origin'")

    origin4 = as_float_list(header4d["space origin"])
    if len(origin4) < 3:
        raise ValueError(f"Invalid 4D origin: {origin4}")
    header3d["space origin"] = [origin4[0], origin4[1], origin4[2]]

    if "kinds" in header4d:
        header3d["kinds"] = list(header4d["kinds"])[:3]

    header3d["endian"] = header4d.get("endian", "little")
    header3d["encoding"] = header4d.get("encoding", "gzip")

    return header3d


def get_4d_header_from_3d_header(header3d, target_path, slice_metadata):
    if "space directions" not in header3d:
        raise ValueError("Input 3D NRRD is missing 'space directions'")
    if "space origin" not in header3d:
        raise ValueError("Input 3D NRRD is missing 'space origin'")

    header4d = {}
    header4d["type"] = header3d.get("type", "float")
    header4d["dimension"] = 4
    header4d["space dimension"] = 4

    dirs3 = header3d["space directions"]
    dirs4 = []
    for d in range(3):
        v = as_float_list(dirs3[d])
        if len(v) < 3:
            raise ValueError(f"Invalid 3D direction vector: {v}")
        dirs4.append([v[0], v[1], v[2], 0.0])
    dirs4.append([0.0, 0.0, 0.0, 1.0])
    header4d["space directions"] = dirs4

    origin3 = as_float_list(header3d["space origin"])
    if len(origin3) < 3:
        raise ValueError(f"Invalid 3D origin: {origin3}")
    header4d["space origin"] = [origin3[0], origin3[1], origin3[2], 0.0]

    header4d["kinds"] = ["domain", "domain", "domain", "domain"]
    header4d["endian"] = header3d.get("endian", "little")
    header4d["encoding"] = "gzip"

    for i, item in enumerate(slice_metadata):
        header4d[f"slice{i}_pixelInterpretation"] = item.get(
            "pixelInterpretation", "SignedDistanceFunction"
        )
        header4d[f"slice{i}_region"] = str(item["region"])
        header4d[f"slice{i}_targetPath"] = item.get("targetPath", target_path)

        site = item.get("site", None)
        if site is not None and str(site) != "":
            header4d[f"slice{i}_site"] = str(site)

    return header4d


def validate_same_geometry(reference_data, reference_header, data, header, path):
    if data.shape != reference_data.shape:
        raise ValueError(
            f"Shape mismatch for {path}: expected {reference_data.shape}, got {data.shape}"
        )

    for key in ["space directions", "space origin"]:
        if key not in reference_header or key not in header:
            raise ValueError(f"Missing '{key}' in {path}")

        a = np.asarray(reference_header[key], dtype=float)
        b = np.asarray(header[key], dtype=float)

        if a.shape != b.shape or not np.allclose(a, b, rtol=0.0, atol=1e-6):
            raise ValueError(
                f"Geometry mismatch in {path}: {key} differs\n"
                f"reference={a}\ncurrent={b}"
            )


def load_slice_metadata_json(path):
    """
    Optional write-mode metadata file.

    Accepted formats:
      1. List of objects in output slice order:
         [
           {"region": "CALC", "site": "ProximalLAD"},
           {"region": "LRNC", "site": "MidLAD"}
         ]

      2. Object with "slices":
         {"slices": [...]}

      3. Mapping from file stem/site/region to metadata:
         {
           "ProximalLAD": {"region": "CALC", "site": "ProximalLAD"},
           "MidLAD": {"region": "LRNC", "site": "MidLAD"}
         }

    Without this file, writer preserves existing region-only behavior unless
    --sites is provided.
    """
    if path is None:
        return None

    with Path(path).open("r", encoding="utf-8") as f:
        payload = json.load(f)

    if isinstance(payload, dict) and "slices" in payload:
        payload = payload["slices"]

    if isinstance(payload, list):
        out = []
        for i, item in enumerate(payload):
            if not isinstance(item, dict):
                raise ValueError(
                    "List-form slice metadata must contain objects/dicts only."
                )
            out.append(
                {
                    "region": str(item.get("region", item.get("name", f"slice{i}"))),
                    "site": None if item.get("site") is None else str(item.get("site")),
                    "pixelInterpretation": str(
                        item.get("pixelInterpretation", "SignedDistanceFunction")
                    ),
                    "targetPath": str(item.get("targetPath", "")),
                }
            )
        return out

    if isinstance(payload, dict):
        out = []
        for key, value in payload.items():
            if isinstance(value, dict):
                out.append(
                    {
                        "region": str(value.get("region", key)),
                        "site": None if value.get("site") is None else str(value.get("site")),
                        "pixelInterpretation": str(
                            value.get("pixelInterpretation", "SignedDistanceFunction")
                        ),
                        "targetPath": str(value.get("targetPath", "")),
                    }
                )
            else:
                out.append(
                    {
                        "region": str(value),
                        "site": str(key),
                        "pixelInterpretation": "SignedDistanceFunction",
                        "targetPath": "",
                    }
                )
        return out

    raise ValueError(f"Unsupported slice metadata JSON format: {path}")


def infer_slice_metadata_from_input_dir(regions, site_names=None):
    """
    Preserve existing behavior when site_names is not provided:
      - require files named by region in REGIONS order.

    New site-aware behavior:
      - if site_names is provided, use files named <site>.nrrd.
      - output order follows the provided site_names list. This makes no
        assumption about lexical order of site names.
    """
    if site_names:
        metadata = []
        for i, site in enumerate(site_names):
            site = str(site)
            region = regions[i] if i < len(regions) else site
            metadata.append(
                {
                    "region": str(region),
                    "site": site,
                    "pixelInterpretation": "SignedDistanceFunction",
                    "targetPath": "",
                }
            )
        return metadata

    return [
        {
            "region": str(region),
            "site": None,
            "pixelInterpretation": "SignedDistanceFunction",
            "targetPath": "",
        }
        for region in regions
    ]


def find_input_nrrd_for_slice(input_dir, item):
    """
    First-level check:
      - If site is present, use site name as file name.
      - Otherwise use region name as before.
    """
    input_dir = Path(input_dir)
    site = item.get("site", None)

    if site is not None and str(site) != "":
        path = input_dir / f"{site}.nrrd"
        if path.exists():
            return path
        raise FileNotFoundError(
            f"Missing required site-named file: {path}. "
            f"When slice site is used, input files must be named by site."
        )

    region = item["region"]
    path = input_dir / f"{region}.nrrd"
    if path.exists():
        return path

    raise FileNotFoundError(f"Missing required file: {path}")


class MultiNrrdReader:
    def __init__(self, input_path):
        self.input_path = Path(input_path)

    def split(self, output_dir):
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        data, header = nrrd.read(str(self.input_path))
        if data.ndim != 4:
            raise ValueError(f"Expected 4D multi-NRRD, got shape {data.shape}")

        n_slices = data.shape[3]
        metadata = parse_slice_metadata(header, n_slices)
        header3d = get_3d_header_from_4d_header(header)

        print(f"Input: {self.input_path}")
        print(f"Shape: {data.shape}")
        print(f"Number of blocks: {n_slices}")
        print(f"3D output header space directions: {header3d['space directions']}")
        print(f"3D output header space origin: {header3d['space origin']}")

        for i in range(n_slices):
            item = metadata[i]
            region = item["region"]
            site = item.get("site", None)
            primary_name = item["primaryName"]

            filename = sanitize_filename(primary_name) + ".nrrd"
            output_path = output_dir / filename

            block = np.ascontiguousarray(data[:, :, :, i], dtype=np.float32)
            nrrd.write(str(output_path), block, header=header3d)

            if site is not None:
                print(
                    f"Wrote slice {i} site={site} region={region}: {output_path}"
                )
            else:
                print(f"Wrote slice {i} region={region}: {output_path}")

        return output_dir


class MultiNrrdWriter:
    def __init__(
        self,
        input_dir,
        target_path,
        regions=None,
        site_names=None,
        slice_metadata_json=None,
    ):
        self.input_dir = Path(input_dir)
        self.target_path = target_path
        self.regions = list(regions or REGIONS)
        self.site_names = None if site_names is None else list(site_names)
        self.slice_metadata = load_slice_metadata_json(slice_metadata_json)

        if self.slice_metadata is None:
            self.slice_metadata = infer_slice_metadata_from_input_dir(
                regions=self.regions,
                site_names=self.site_names,
            )

        for item in self.slice_metadata:
            if not item.get("targetPath"):
                item["targetPath"] = self.target_path
            if not item.get("pixelInterpretation"):
                item["pixelInterpretation"] = "SignedDistanceFunction"
            if "region" not in item or item["region"] is None:
                item["region"] = item.get("site", "slice")

    def write(self, output_dir, output_name="composition.multi.nrrd"):
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        output_path = output_dir / output_name

        arrays = []
        reference_data = None
        reference_header = None

        for item in self.slice_metadata:
            path = find_input_nrrd_for_slice(self.input_dir, item)

            data, header = nrrd.read(str(path))
            if data.ndim != 3:
                raise ValueError(f"Expected 3D NRRD for {path}, got shape {data.shape}")

            if reference_data is None:
                reference_data = data
                reference_header = header
            else:
                validate_same_geometry(reference_data, reference_header, data, header, path)

            arrays.append(np.asarray(data, dtype=np.float32))

            site = item.get("site", None)
            region = item["region"]
            if site is not None and str(site) != "":
                print(f"Read site={site} region={region} from {path}: shape={data.shape}")
            else:
                print(f"Read region={region} from {path}: shape={data.shape}")

        stacked = np.stack(arrays, axis=3)
        stacked = np.ascontiguousarray(stacked.astype(np.float32))

        header4d = get_4d_header_from_3d_header(
            reference_header,
            target_path=self.target_path,
            slice_metadata=self.slice_metadata,
        )

        nrrd.write(str(output_path), stacked, header=header4d)
        print(f"Wrote composition multi-NRRD: {output_path}")
        print(f"Shape: {stacked.shape}")

        return output_path


def parse_args(argv):
    parser = argparse.ArgumentParser(
        description="Read/write composition .multi.nrrd files."
    )
    subparsers = parser.add_subparsers(dest="mode", required=True)

    read_parser = subparsers.add_parser("read")
    read_parser.add_argument("--input", required=True, help="Input .multi.nrrd")
    read_parser.add_argument("--output-dir", required=True, help="Output directory")

    write_parser = subparsers.add_parser("write")
    write_parser.add_argument("--input-dir", required=True, help="Directory with region/site NRRDs")
    write_parser.add_argument("--output-dir", required=True, help="Output directory")
    write_parser.add_argument("--target-path", required=True, help="e.g. LeftCoronary")
    write_parser.add_argument(
        "--output-name",
        default="composition.multi.nrrd",
        help="Default: composition.multi.nrrd",
    )
    write_parser.add_argument(
        "--sites",
        nargs="+",
        default=None,
        help=(
            "Optional site names to write as slice{i}_site. "
            "When provided, input files are read as <site>.nrrd in this order."
        ),
    )
    write_parser.add_argument(
        "--slice-metadata-json",
        default=None,
        help=(
            "Optional JSON file describing slices. Supports region, site, "
            "pixelInterpretation, and targetPath. If site is present, input "
            "file name is <site>.nrrd."
        ),
    )

    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv or sys.argv[1:])

    if args.mode == "read":
        MultiNrrdReader(args.input).split(args.output_dir)
    elif args.mode == "write":
        MultiNrrdWriter(
            args.input_dir,
            args.target_path,
            site_names=args.sites,
            slice_metadata_json=args.slice_metadata_json,
        ).write(
            args.output_dir,
            args.output_name,
        )


if __name__ == "__main__":
    main()
