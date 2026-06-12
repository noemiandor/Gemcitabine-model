#!/usr/bin/env python3
"""Generate a six-timepoint class-contour overlay panel.

This wrapper calls scripts/generate_class_colored_overlay.py for each requested
elapsed time, then assembles the generated overlays into one labeled image.
"""

from __future__ import annotations

import argparse
import math
import subprocess
import sys
from pathlib import Path
from typing import Sequence

from PIL import Image, ImageDraw, ImageFont


DEFAULT_HOURS = (2, 12, 24, 36, 48, 60)


def parse_args() -> argparse.Namespace:
    script_dir = Path(__file__).resolve().parent
    project_root = script_dir.parent
    analysis_dir = project_root / "Final_Tracking_analysis"

    parser = argparse.ArgumentParser(
        description="Generate D5_2 class-contour overlays and a timecourse panel."
    )
    parser.add_argument("--well", default="D5_2")
    parser.add_argument(
        "--hours",
        type=int,
        nargs="+",
        default=list(DEFAULT_HOURS),
        help="Elapsed hours to render. Each must be divisible by 2.",
    )
    parser.add_argument(
        "--analysis-dir",
        default=str(analysis_dir),
        help="Final_Tracking_analysis directory.",
    )
    parser.add_argument(
        "--overlay-script",
        default=str(script_dir / "generate_class_colored_overlay.py"),
        help="Path to the single-frame overlay script.",
    )
    parser.add_argument(
        "--output-dir",
        default=None,
        help=(
            "Directory for individual overlays and the assembled panel. Defaults "
            "to Final_Tracking_analysis/Overlays/<well>_timecourse."
        ),
    )
    parser.add_argument(
        "--classification-subset",
        default=None,
        help="Optional classifier CSV subset to pass through for faster rendering.",
    )
    parser.add_argument(
        "--max-match-distance",
        type=float,
        default=15.0,
        help="Pixel threshold passed to the single-frame overlay script.",
    )
    parser.add_argument(
        "--panel-width",
        type=int,
        default=900,
        help="Width of each panel in the assembled timecourse image.",
    )
    parser.add_argument(
        "--columns",
        type=int,
        default=3,
        help="Number of columns in the assembled timecourse image.",
    )
    parser.add_argument(
        "--crop-fraction",
        type=float,
        default=1.0,
        help=(
            "Fraction of the full field-of-view area to show in the assembled "
            "panel. The crop is centered and keeps the original aspect ratio."
        ),
    )
    parser.add_argument(
        "--crop-center-x",
        type=float,
        default=0.5,
        help="Horizontal crop center as a fraction of image width.",
    )
    parser.add_argument(
        "--crop-center-y",
        type=float,
        default=0.5,
        help="Vertical crop center as a fraction of image height.",
    )
    parser.add_argument(
        "--skip-render",
        action="store_true",
        help="Assemble existing overlay PNGs without rerunning the overlay script.",
    )
    return parser.parse_args()


def elapsed_label(hours: int) -> str:
    days, hour = divmod(hours, 24)
    return f"{days:02d}d{hour:02d}h00m"


def frame_for_hours(hours: int) -> int:
    if hours < 0 or hours % 2 != 0:
        raise ValueError(f"Elapsed time must be a nonnegative multiple of 2 h: {hours}")
    return hours // 2 + 1


def overlay_filename(well: str, hours: int, frame: int) -> str:
    return f"{well}_{hours:03d}h_t{frame:02d}_class_contour_overlay.png"


def run_overlay_script(
    overlay_script: Path,
    analysis_dir: Path,
    well: str,
    hours: int,
    output_path: Path,
    classification_subset: str | None,
    max_match_distance: float,
) -> None:
    frame = frame_for_hours(hours)
    command = [
        sys.executable,
        str(overlay_script),
        "--analysis-dir",
        str(analysis_dir),
        "--well",
        well,
        "--frame",
        str(frame),
        "--classification-time",
        elapsed_label(hours),
        "--max-match-distance",
        str(max_match_distance),
        "--output",
        str(output_path),
    ]
    if classification_subset:
        command.extend(["--classification-subset", classification_subset])

    subprocess.run(command, check=True)


def load_font(size: int) -> ImageFont.ImageFont:
    for font_path in (
        "/System/Library/Fonts/Supplemental/Arial.ttf",
        "/Library/Fonts/Arial.ttf",
        "/System/Library/Fonts/Helvetica.ttc",
        "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",
    ):
        try:
            return ImageFont.truetype(font_path, size=size)
        except OSError:
            continue
    return ImageFont.load_default()


def resize_panel(image: Image.Image, width: int) -> Image.Image:
    ratio = width / image.width
    height = int(round(image.height * ratio))
    return image.resize((width, height), Image.Resampling.LANCZOS)


def crop_field_of_view(
    image: Image.Image,
    crop_fraction: float,
    center_x_fraction: float,
    center_y_fraction: float,
) -> Image.Image:
    if not 0 < crop_fraction <= 1:
        raise ValueError("--crop-fraction must be > 0 and <= 1")
    if not 0 <= center_x_fraction <= 1:
        raise ValueError("--crop-center-x must be between 0 and 1")
    if not 0 <= center_y_fraction <= 1:
        raise ValueError("--crop-center-y must be between 0 and 1")
    if crop_fraction == 1:
        return image

    side_scale = math.sqrt(crop_fraction)
    crop_width = max(1, int(round(image.width * side_scale)))
    crop_height = max(1, int(round(image.height * side_scale)))
    center_x = image.width * center_x_fraction
    center_y = image.height * center_y_fraction
    left = int(round(center_x - crop_width / 2))
    top = int(round(center_y - crop_height / 2))
    left = min(max(left, 0), image.width - crop_width)
    top = min(max(top, 0), image.height - crop_height)

    return image.crop((left, top, left + crop_width, top + crop_height))


def assemble_panel(
    overlay_paths: Sequence[Path],
    hours: Sequence[int],
    well: str,
    output_path: Path,
    panel_width: int,
    columns: int,
    crop_fraction: float,
    crop_center_x: float,
    crop_center_y: float,
) -> None:
    if columns <= 0:
        raise ValueError("--columns must be positive")
    if len(overlay_paths) != len(hours):
        raise ValueError("overlay_paths and hours must have the same length")

    panels = [
        resize_panel(
            crop_field_of_view(
                Image.open(path).convert("RGB"),
                crop_fraction,
                crop_center_x,
                crop_center_y,
            ),
            panel_width,
        )
        for path in overlay_paths
    ]
    panel_height = panels[0].height
    label_height = max(48, panel_width // 16)
    gutter = max(16, panel_width // 40)
    rows = math.ceil(len(panels) / columns)

    canvas_width = columns * panel_width + (columns + 1) * gutter
    canvas_height = rows * (panel_height + label_height) + (rows + 1) * gutter
    canvas = Image.new("RGB", (canvas_width, canvas_height), "white")
    draw = ImageDraw.Draw(canvas)
    font = load_font(max(22, panel_width // 30))
    small_font = load_font(max(16, panel_width // 45))

    for index, (panel, hour) in enumerate(zip(panels, hours)):
        row, col = divmod(index, columns)
        x = gutter + col * (panel_width + gutter)
        y = gutter + row * (panel_height + label_height + gutter)
        title = f"{hour} h"
        subtitle = f"{well} frame {frame_for_hours(hour)}"

        title_box = draw.textbbox((0, 0), title, font=font)
        subtitle_box = draw.textbbox((0, 0), subtitle, font=small_font)
        draw.text(
            (x + (panel_width - (title_box[2] - title_box[0])) / 2, y),
            title,
            fill="black",
            font=font,
        )
        draw.text(
            (
                x + (panel_width - (subtitle_box[2] - subtitle_box[0])) / 2,
                y + max(26, panel_width // 36),
            ),
            subtitle,
            fill=(70, 70, 70),
            font=small_font,
        )
        canvas.paste(panel, (x, y + label_height))

    output_path.parent.mkdir(parents=True, exist_ok=True)
    canvas.save(output_path)


def main() -> None:
    args = parse_args()
    analysis_dir = Path(args.analysis_dir)
    overlay_script = Path(args.overlay_script)
    output_dir = (
        Path(args.output_dir)
        if args.output_dir is not None
        else analysis_dir / "Overlays" / f"{args.well}_timecourse"
    )
    output_dir.mkdir(parents=True, exist_ok=True)

    hours = list(args.hours)
    overlay_paths = []
    for hour in hours:
        frame = frame_for_hours(hour)
        output_path = output_dir / overlay_filename(args.well, hour, frame)
        overlay_paths.append(output_path)
        if not args.skip_render:
            print(
                f"rendering {args.well} {hour} h -> frame {frame}: {output_path}",
                flush=True,
            )
            run_overlay_script(
                overlay_script=overlay_script,
                analysis_dir=analysis_dir,
                well=args.well,
                hours=hour,
                output_path=output_path,
                classification_subset=args.classification_subset,
                max_match_distance=args.max_match_distance,
            )

    panel_path = output_dir / (
        f"{args.well}_timecourse_"
        + "_".join(f"{hour:03d}h" for hour in hours)
        + f"_crop{args.crop_fraction:g}.png"
    )
    assemble_panel(
        overlay_paths=overlay_paths,
        hours=hours,
        well=args.well,
        output_path=panel_path,
        panel_width=args.panel_width,
        columns=args.columns,
        crop_fraction=args.crop_fraction,
        crop_center_x=args.crop_center_x,
        crop_center_y=args.crop_center_y,
    )

    print(f"wrote panel: {panel_path}")
    print("wrote overlays:")
    for path in overlay_paths:
        print(f"  {path}")


if __name__ == "__main__":
    main()
