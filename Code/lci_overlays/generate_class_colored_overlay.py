#!/usr/bin/env python3
"""Generate a class-colored mask contour overlay on a raw live-cell TIFF frame.

This script assumes the layout used by:
Final_Tracking_analysis/
  Images_40Frames/
  Masks/
  Cell_classification_CSVs/

Mask frame numbers are treated as 1-based: frame 24 uses stack index 23.
The classification CSV contains 2-hour timepoint image names. By default this
script maps mask/stack frame N to elapsed time (N - 1) * 2 hours. Use
--classification-time to override that mapping.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import Iterable, Tuple

import numpy as np
import tifffile as tiff
from PIL import Image
from scipy import ndimage as ndi
from scipy.spatial import cKDTree


DEFAULT_ANALYSIS_DIR = "Final_Tracking_analysis"

PALETTE = {
    "Alive": np.array([34, 197, 94], dtype=np.uint8),
    "Dead": np.array([239, 68, 68], dtype=np.uint8),
    "Transitional": np.array([245, 158, 11], dtype=np.uint8),
    "UnStained": np.array([59, 130, 246], dtype=np.uint8),
    "Unknown": np.array([168, 168, 168], dtype=np.uint8),
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Render a class-colored segmentation overlay for one well/frame."
    )
    parser.add_argument("--analysis-dir", default=DEFAULT_ANALYSIS_DIR)
    parser.add_argument("--well", default="D5_2")
    parser.add_argument(
        "--frame",
        type=int,
        default=24,
        help="1-based mask/stack frame number, e.g. 24 for D5_2_t24.tiff.",
    )
    parser.add_argument(
        "--classification-time",
        default=None,
        help="Classifier image time label, e.g. 05d18h00m. Defaults from frame.",
    )
    parser.add_argument(
        "--frame-interval-hours",
        type=int,
        default=2,
        help="Hours between frames in Images_40Frames, used for classifier mapping.",
    )
    parser.add_argument("--output", default=None)
    parser.add_argument(
        "--fill-weight",
        type=float,
        default=0.0,
        help=argparse.SUPPRESS,
    )
    parser.add_argument(
        "--classification-subset",
        default=None,
        help=(
            "Optional CSV subset for the well. If omitted, the row classifier CSV "
            "is scanned directly."
        ),
    )
    parser.add_argument(
        "--max-match-distance",
        type=float,
        default=15.0,
        help=(
            "Maximum pixel distance between a mask-component centroid and classifier "
            "object box center before the contour is marked Unknown."
        ),
    )
    return parser.parse_args()


def elapsed_label(hours: int) -> str:
    days, hour = divmod(hours, 24)
    return f"{days:02d}d{hour:02d}h00m"


def read_raw_frame(stack_path: Path, frame_number: int) -> np.ndarray:
    stack = tiff.imread(stack_path)
    frame_index = frame_number - 1
    if frame_index < 0:
        raise ValueError("--frame must be 1 or greater")

    if stack.ndim == 4:
        raw = stack[frame_index]
    elif stack.ndim == 3 and stack.shape[-1] in (3, 4):
        # Single RGB/RGBA image.
        if frame_number != 1:
            raise ValueError(f"{stack_path} contains one RGB image, not a stack")
        raw = stack[..., :3]
    elif stack.ndim == 3:
        raw = stack[frame_index]
        raw = np.stack([raw] * 3, axis=-1)
    elif stack.ndim == 2:
        if frame_number != 1:
            raise ValueError(f"{stack_path} contains one grayscale image, not a stack")
        raw = np.stack([stack] * 3, axis=-1)
    else:
        raise ValueError(f"Unsupported TIFF shape for {stack_path}: {stack.shape}")

    if raw.ndim == 2:
        raw = np.stack([raw] * 3, axis=-1)
    if raw.shape[-1] == 4:
        raw = raw[..., :3]

    if raw.dtype == np.uint8:
        return raw

    raw_float = raw.astype(np.float32)
    lo, hi = np.percentile(raw_float, [0.1, 99.9])
    if hi <= lo:
        lo, hi = float(raw_float.min()), float(raw_float.max())
    return np.clip(255 * (raw_float - lo) / (hi - lo + 1e-9), 0, 255).astype(np.uint8)


def read_binary_mask(mask_path: Path) -> np.ndarray:
    return np.array(Image.open(mask_path).convert("L")) > 0


def iter_filtered_classifier_rows(
    csv_path: Path,
    well: str,
    classification_time: str,
) -> Iterable[list[str]]:
    well_bytes = well.encode()
    time_bytes = classification_time.encode()
    with csv_path.open("rb") as handle:
        header_line = handle.readline().decode("utf-8-sig", errors="replace")
        header = next(csv.reader([header_line]))
        for raw_line in handle:
            if well_bytes not in raw_line or time_bytes not in raw_line:
                continue
            line = raw_line.decode("utf-8-sig", errors="replace")
            row = next(csv.reader([line]))
            if len(row) == len(header):
                yield row


def iter_subset_classifier_rows(
    subset_path: Path,
    classification_time: str,
) -> Iterable[list[str]]:
    with subset_path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.reader(handle)
        for row in reader:
            if len(row) >= 5 and classification_time in row[0]:
                yield row


def read_classifier_objects(
    classifier_csv: Path,
    well: str,
    classification_time: str,
    subset_path: Path | None = None,
) -> Tuple[np.ndarray, list[str]]:
    rows: Iterable[list[str]]
    if subset_path is not None:
        rows = iter_subset_classifier_rows(subset_path, classification_time)
    else:
        rows = iter_filtered_classifier_rows(classifier_csv, well, classification_time)

    centers = []
    classes = []
    for row in rows:
        xmin, xmax, ymin, ymax = (float(value) for value in row[5:9])
        centers.append(((xmin + xmax) / 2.0, (ymin + ymax) / 2.0))
        classes.append(row[4].strip())
    if not centers:
        source = subset_path if subset_path is not None else classifier_csv
        raise ValueError(
            f"No classifier rows for {well} at {classification_time} in {source}"
        )
    return np.asarray(centers, dtype=float), classes


def edge_mask(mask: np.ndarray) -> np.ndarray:
    edge = np.zeros_like(mask, dtype=bool)
    edge[1:, :] |= mask[1:, :] != mask[:-1, :]
    edge[:-1, :] |= mask[:-1, :] != mask[1:, :]
    edge[:, 1:] |= mask[:, 1:] != mask[:, :-1]
    edge[:, :-1] |= mask[:, :-1] != mask[:, 1:]
    return edge


def render_overlay(
    raw: np.ndarray,
    mask: np.ndarray,
    classifier_centers: np.ndarray,
    classifier_classes: list[str],
    max_match_distance: float,
) -> Tuple[np.ndarray, int, int]:
    labels, component_count = ndi.label(mask)
    tree = cKDTree(classifier_centers)

    output = np.zeros((raw.shape[0], raw.shape[1], 4), dtype=np.uint8)
    output[..., :3] = raw
    output[..., 3] = 255

    matched_components = 0

    for component_id in range(1, component_count + 1):
        component_mask = labels == component_id
        if not component_mask.any():
            continue
        center_y, center_x = ndi.center_of_mass(component_mask)
        distance, classifier_index = tree.query([center_x, center_y], k=1)
        if distance <= max_match_distance:
            phenotype = classifier_classes[int(classifier_index)]
            matched_components += 1
        else:
            phenotype = "Unknown"
        color = PALETTE.get(phenotype, PALETTE["Unknown"])

        component_edge = edge_mask(component_mask)
        output[component_edge] = [int(color[0]), int(color[1]), int(color[2]), 255]

    return output, component_count, matched_components


def main() -> None:
    args = parse_args()
    analysis_dir = Path(args.analysis_dir)
    well = args.well
    row = well[0]
    frame = args.frame

    classification_time = args.classification_time
    if classification_time is None:
        classification_time = elapsed_label((frame - 1) * args.frame_interval_hours)

    raw_stack = analysis_dir / "Images_40Frames" / f"{well}.tiff"
    mask_path = analysis_dir / "Masks" / well / f"{well}_t{frame}.tiff"
    classifier_csv = (
        analysis_dir / "Cell_classification_CSVs" / f"{row}_Row_Objet_Data.csv"
    )
    output = (
        Path(args.output)
        if args.output is not None
        else Path(f"{well}_t{frame}_class_overlay.png")
    )
    subset_path = Path(args.classification_subset) if args.classification_subset else None

    raw = read_raw_frame(raw_stack, frame)
    mask = read_binary_mask(mask_path)
    classifier_centers, classifier_classes = read_classifier_objects(
        classifier_csv, well, classification_time, subset_path
    )

    overlay, component_count, matched_components = render_overlay(
        raw,
        mask,
        classifier_centers,
        classifier_classes,
        args.max_match_distance,
    )

    output.parent.mkdir(parents=True, exist_ok=True)
    Image.fromarray(overlay, "RGBA").save(output)

    print(f"wrote: {output}")
    print(f"well: {well}")
    print(f"frame: {frame}")
    print(f"classification_time: {classification_time}")
    print(f"mask_components: {component_count}")
    print(f"classifier_objects: {len(classifier_classes)}")
    print(f"max_match_distance_px: {args.max_match_distance}")
    print(f"matched_components: {matched_components}")
    print(
        "contour palette: Alive=green Dead=red Transitional=amber "
        "UnStained=blue Unknown=gray"
    )


if __name__ == "__main__":
    main()
