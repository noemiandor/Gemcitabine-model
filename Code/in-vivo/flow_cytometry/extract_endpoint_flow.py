#!/usr/bin/env python3
"""Extract reviewed endpoint DNA-content gates from a FlowJo workspace.

The extractor intentionally reads only FCS metadata. Population membership and
counts are taken from the reviewed FlowJo workspace, whose gates can depend on
FlowJo transformations that should not be approximated in a dependency-light
reader.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import math
import os
import re
import sys
import tempfile
import xml.etree.ElementTree as ET
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from statistics import mean
from typing import Iterable, Mapping, Sequence
from urllib.parse import unquote, urlparse


REQUIRED_CROSSWALK_COLUMNS = (
    "mouse_id",
    "injected_origin",
    "dose_mg_kg",
    "fcs_file",
    "wsp_sample_name",
    "peak_gate_name",
    "reviewed_peak_annotation_n",
    "acquisition_date",
)
REQUIRED_GATE_NAMES = ("2N", "4N", r"mCh+\2N", r"mCh+\4N")
EXPECTED_MANIFEST_COLUMNS = ("input_role", "relative_path", "size_bytes", "sha256")
SENSITIVITY_COLUMNS = (
    "mouse_id",
    "selected_wsp_sample_name",
    "selected_peak_gate_name",
    "selected_peak_annotation_n",
    "paired_wsp_sample_name",
    "paired_peak_gate_name",
    "paired_peak_annotation_n",
    "relationship_status",
)


class ExtractionError(ValueError):
    """Raised when an input violates the reviewed extraction contract."""


@dataclass(frozen=True)
class CrosswalkRow:
    mouse_id: str
    injected_origin: str
    dose_mg_kg: str
    fcs_file: str
    wsp_sample_name: str
    peak_gate_name: str
    reviewed_peak_annotation_n: str
    acquisition_date: str


@dataclass(frozen=True)
class FCSMetadata:
    version: str
    keywords: Mapping[str, str]


@dataclass(frozen=True)
class WorkspaceSample:
    sample_id: str
    name: str
    data_set_basename: str
    total_count: int
    human_cell_enrichment_count: int
    human_count: int
    populations: Mapping[str, int]


@dataclass(frozen=True)
class SensitivityRow:
    mouse_id: str
    selected_wsp_sample_name: str
    selected_peak_gate_name: str
    selected_peak_annotation_n: str
    paired_wsp_sample_name: str
    paired_peak_gate_name: str
    paired_peak_annotation_n: str
    relationship_status: str


PER_SAMPLE_FIELDS = (
    "mouse_id",
    "injected_origin",
    "dose_mg_kg",
    "fcs_file",
    "wsp_sample_name",
    "wsp_sample_id",
    "fcs_total_count",
    "wsp_total_count",
    "human_cell_enrichment_count",
    "human_cell_enrichment_pct_total",
    "human_cells_count",
    "human_cells_pct_enrichment",
    "human_cells_pct_total",
    "human_cells_qc_pass",
    "human_cells_qc_note",
    "peak_gate_name",
    "reviewed_peak_annotation_n",
    "peak_count",
    "peak_pct_human_cells",
    "two_n_count",
    "two_n_pct_human_cells",
    "four_n_count",
    "four_n_pct_human_cells",
    "mch_pos_two_n_count",
    "mch_pos_two_n_pct_human_cells",
    "mch_pos_four_n_count",
    "mch_pos_four_n_pct_human_cells",
    "summed_mch_pos_two_n_four_n_gate_count",
    "mch_pos_two_n_pct_of_summed_two_n_four_n_gate_counts",
    "mch_pos_four_n_pct_of_summed_two_n_four_n_gate_counts",
    "fcs_version",
    "fcs_cytometer",
    "dna_content_parameter_index",
    "dna_content_parameter_name",
    "crosswalk_acquisition_date",
    "fcs_date",
    "fcs_begin_time",
    "fcs_end_time",
)


GROUP_FIELDS = (
    "group_level",
    "injected_origin",
    "dose_mg_kg",
    "n_samples",
    "n_samples_human_cells_qc_pass",
    "total_human_cells_count",
    "total_summed_mch_pos_two_n_four_n_gate_count",
    "reviewed_peak_gate_labels",
    "mean_peak_pct_human_cells_all_samples",
    "mean_peak_pct_human_cells_qc_pass",
    "mean_two_n_pct_human_cells",
    "mean_four_n_pct_human_cells",
    "mean_mch_pos_four_n_pct_of_summed_two_n_four_n_gate_counts",
    "pooled_mch_pos_four_n_pct_of_summed_two_n_four_n_gate_counts",
)


SENSITIVITY_OUTPUT_FIELDS = (
    "mouse_id",
    "selected_wsp_sample_name",
    "selected_peak_gate_name",
    "selected_peak_annotation_n",
    "selected_human_cells_count",
    "selected_peak_count",
    "selected_peak_pct_human_cells",
    "paired_wsp_sample_name",
    "paired_peak_gate_name",
    "paired_peak_annotation_n",
    "paired_human_cells_count",
    "paired_peak_count",
    "paired_peak_pct_human_cells",
    "relationship_status",
    "paired_fcs_imported",
)


def _clean_field(value: object, field: str, line_number: int) -> str:
    text = "" if value is None else str(value).strip()
    if not text:
        raise ExtractionError(
            f"crosswalk line {line_number}: required field {field!r} is empty"
        )
    return text


def read_crosswalk(path: Path, expected_samples: int) -> list[CrosswalkRow]:
    try:
        with path.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            fields = tuple(reader.fieldnames or ())
            missing = [name for name in REQUIRED_CROSSWALK_COLUMNS if name not in fields]
            if missing:
                raise ExtractionError(
                    "crosswalk is missing required column(s): " + ", ".join(missing)
                )
            rows = []
            for line_number, raw in enumerate(reader, start=2):
                values = {
                    name: _clean_field(raw.get(name), name, line_number)
                    for name in REQUIRED_CROSSWALK_COLUMNS
                }
                if values["injected_origin"] not in {"2N", "4N"}:
                    raise ExtractionError(
                        f"crosswalk line {line_number}: injected_origin must be 2N or 4N"
                    )
                try:
                    dose = float(values["dose_mg_kg"])
                except ValueError as exc:
                    raise ExtractionError(
                        f"crosswalk line {line_number}: dose_mg_kg is not numeric"
                    ) from exc
                if not math.isfinite(dose) or dose < 0:
                    raise ExtractionError(
                        f"crosswalk line {line_number}: dose_mg_kg must be finite and nonnegative"
                    )
                if Path(values["fcs_file"]).is_absolute():
                    raise ExtractionError(
                        f"crosswalk line {line_number}: fcs_file must be relative to --fcs-dir"
                    )
                if ".." in Path(values["fcs_file"]).parts:
                    raise ExtractionError(
                        f"crosswalk line {line_number}: fcs_file cannot traverse outside --fcs-dir"
                    )
                try:
                    reviewed_peak = float(values["reviewed_peak_annotation_n"])
                except ValueError as exc:
                    raise ExtractionError(
                        f"crosswalk line {line_number}: reviewed_peak_annotation_n is not numeric"
                    ) from exc
                if not math.isfinite(reviewed_peak) or reviewed_peak <= 0:
                    raise ExtractionError(
                        f"crosswalk line {line_number}: reviewed_peak_annotation_n must be "
                        "finite and positive"
                    )
                try:
                    datetime.strptime(values["acquisition_date"], "%Y-%m-%d")
                except ValueError as exc:
                    raise ExtractionError(
                        f"crosswalk line {line_number}: acquisition_date must use YYYY-MM-DD"
                    ) from exc
                rows.append(CrosswalkRow(**values))
    except OSError as exc:
        raise ExtractionError(f"cannot read crosswalk {path}: {exc}") from exc

    if len(rows) != expected_samples:
        raise ExtractionError(
            f"crosswalk has {len(rows)} rows; expected exactly {expected_samples}"
        )
    for field in ("mouse_id", "fcs_file", "wsp_sample_name"):
        values = [getattr(row, field) for row in rows]
        duplicates = sorted({value for value in values if values.count(value) > 1})
        if duplicates:
            raise ExtractionError(
                f"crosswalk {field} values must be unique; duplicate(s): "
                + ", ".join(duplicates)
            )
    return rows


def read_sensitivity_table(path: Path) -> list[SensitivityRow]:
    try:
        with path.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if tuple(reader.fieldnames or ()) != SENSITIVITY_COLUMNS:
                raise ExtractionError(
                    "paired-sensitivity columns must be exactly: "
                    + ", ".join(SENSITIVITY_COLUMNS)
                )
            rows = []
            for line_number, raw in enumerate(reader, start=2):
                values = {
                    name: _clean_field(raw.get(name), name, line_number)
                    for name in SENSITIVITY_COLUMNS
                }
                for field in (
                    "selected_peak_annotation_n",
                    "paired_peak_annotation_n",
                ):
                    try:
                        value = float(values[field])
                    except ValueError as exc:
                        raise ExtractionError(
                            f"paired-sensitivity line {line_number}: {field} is not numeric"
                        ) from exc
                    if not math.isfinite(value) or value <= 0:
                        raise ExtractionError(
                            f"paired-sensitivity line {line_number}: {field} must be "
                            "finite and positive"
                        )
                if values["relationship_status"] != "unresolved":
                    raise ExtractionError(
                        f"paired-sensitivity line {line_number}: relationship_status "
                        "must remain unresolved until primary metadata establishes it"
                    )
                rows.append(SensitivityRow(**values))
    except OSError as exc:
        raise ExtractionError(f"cannot read paired-sensitivity table {path}: {exc}") from exc
    if not rows:
        raise ExtractionError("paired-sensitivity table must contain at least one row")
    for field in ("mouse_id", "paired_wsp_sample_name"):
        values = [getattr(row, field) for row in rows]
        duplicates = sorted({value for value in values if values.count(value) > 1})
        if duplicates:
            raise ExtractionError(
                f"paired-sensitivity {field} values must be unique; duplicate(s): "
                + ", ".join(duplicates)
            )
    return rows


def _parse_header_offset(raw: bytes, label: str) -> int:
    text = raw.decode("ascii", errors="strict").strip()
    if not text:
        raise ExtractionError(f"FCS header has an empty {label} offset")
    try:
        return int(text)
    except ValueError as exc:
        raise ExtractionError(f"FCS header has an invalid {label} offset: {text!r}") from exc


def _parse_fcs_text(raw: bytes) -> dict[str, str]:
    try:
        text = raw.decode("latin-1")
    except UnicodeDecodeError as exc:  # pragma: no cover - latin-1 accepts all bytes
        raise ExtractionError("FCS TEXT segment is not decodable") from exc
    if len(text) < 3:
        raise ExtractionError("FCS TEXT segment is too short")
    delimiter = text[0]
    if not delimiter or delimiter.isalnum():
        raise ExtractionError("FCS TEXT segment has an invalid delimiter")

    # FCS uses a doubled delimiter to quote the delimiter within a token.
    tokens: list[str] = []
    token: list[str] = []
    index = 1
    while index < len(text):
        char = text[index]
        if char != delimiter:
            token.append(char)
            index += 1
            continue
        if index + 1 < len(text) and text[index + 1] == delimiter:
            token.append(delimiter)
            index += 2
            continue
        tokens.append("".join(token))
        token = []
        index += 1
    if token:
        tokens.append("".join(token))
    if tokens and tokens[-1] == "":
        tokens.pop()
    if len(tokens) % 2:
        raise ExtractionError("FCS TEXT segment does not contain key/value pairs")

    keywords: dict[str, str] = {}
    for index in range(0, len(tokens), 2):
        key, value = tokens[index], tokens[index + 1]
        if not key:
            raise ExtractionError("FCS TEXT segment contains an empty keyword")
        if key in keywords and keywords[key] != value:
            raise ExtractionError(f"FCS TEXT segment repeats keyword {key!r}")
        keywords[key] = value.strip()
    return keywords


def read_fcs_metadata(path: Path) -> FCSMetadata:
    try:
        size = path.stat().st_size
        with path.open("rb") as handle:
            header = handle.read(58)
            if len(header) != 58:
                raise ExtractionError(f"FCS file is shorter than its 58-byte header: {path}")
            try:
                version = header[0:6].decode("ascii").strip()
            except UnicodeDecodeError as exc:
                raise ExtractionError(f"FCS version is not ASCII: {path}") from exc
            if not re.fullmatch(r"FCS[23]\.[0-9]", version):
                raise ExtractionError(f"unsupported FCS version {version!r}: {path}")
            text_start = _parse_header_offset(header[10:18], "TEXT start")
            text_end = _parse_header_offset(header[18:26], "TEXT end")
            if text_start < 58 or text_end < text_start or text_end >= size:
                raise ExtractionError(
                    f"FCS TEXT offsets [{text_start}, {text_end}] are outside {path}"
                )
            handle.seek(text_start)
            raw_text = handle.read(text_end - text_start + 1)
    except OSError as exc:
        raise ExtractionError(f"cannot read FCS file {path}: {exc}") from exc

    keywords = _parse_fcs_text(raw_text)
    for keyword in (
        "$TOT",
        "$PAR",
        "$BYTEORD",
        "$DATATYPE",
        "$FIL",
        "$CYT",
        "$DATE",
    ):
        if not keywords.get(keyword):
            raise ExtractionError(f"FCS file {path} lacks required keyword {keyword}")
    try:
        total = int(keywords["$TOT"])
        parameters = int(keywords["$PAR"])
    except ValueError as exc:
        raise ExtractionError(f"FCS file {path} has noninteger $TOT or $PAR") from exc
    if total < 0 or parameters <= 0:
        raise ExtractionError(f"FCS file {path} has invalid $TOT or $PAR")
    return FCSMetadata(version=version, keywords=keywords)


def _tag_name(element: ET.Element) -> str:
    return element.tag.rsplit("}", 1)[-1]


def _direct_children(element: ET.Element, tag: str) -> list[ET.Element]:
    return [child for child in element if _tag_name(child) == tag]


def _one(elements: Sequence[ET.Element], description: str) -> ET.Element:
    if len(elements) != 1:
        raise ExtractionError(f"expected exactly one {description}; found {len(elements)}")
    return elements[0]


def _population_count(element: ET.Element, context: str) -> int:
    raw = element.get("count")
    try:
        count = int(raw or "")
    except ValueError as exc:
        raise ExtractionError(f"{context} has invalid population count {raw!r}") from exc
    if count < 0:
        raise ExtractionError(f"{context} has negative population count")
    return count


def read_workspace(
    path: Path, selected_names: set[str] | None = None
) -> dict[str, WorkspaceSample]:
    try:
        root = ET.parse(path).getroot()
    except (OSError, ET.ParseError) as exc:
        raise ExtractionError(f"cannot parse FlowJo workspace {path}: {exc}") from exc

    sample_list = _one(_direct_children(root, "SampleList"), "SampleList")
    samples: dict[str, WorkspaceSample] = {}
    all_names: set[str] = set()
    sample_ids: set[str] = set()
    for sample in _direct_children(sample_list, "Sample"):
        data_set = _one(_direct_children(sample, "DataSet"), "DataSet in Sample")
        sample_node = _one(_direct_children(sample, "SampleNode"), "SampleNode in Sample")
        sample_id = (sample_node.get("sampleID") or "").strip()
        data_set_id = (data_set.get("sampleID") or "").strip()
        name = (sample_node.get("name") or "").strip()
        if not sample_id or sample_id != data_set_id:
            raise ExtractionError(
                f"workspace sample {name!r} has inconsistent DataSet/SampleNode IDs"
            )
        if not name:
            raise ExtractionError(f"workspace sampleID {sample_id!r} has no name")
        if name in all_names:
            raise ExtractionError(f"workspace contains duplicate SampleNode name {name!r}")
        all_names.add(name)
        if sample_id in sample_ids:
            raise ExtractionError(f"workspace contains duplicate sampleID {sample_id!r}")
        sample_ids.add(sample_id)

        uri = (data_set.get("uri") or "").strip()
        parsed_uri = urlparse(uri)
        data_set_basename = Path(unquote(parsed_uri.path or uri)).name
        if not data_set_basename:
            raise ExtractionError(f"workspace sample {name!r} has an invalid DataSet URI")

        # A FlowJo workspace often retains exploratory/control samples with
        # incomplete or duplicated gates. The explicit crosswalk defines the
        # reviewed cohort; do not let unrelated nodes enter or block extraction.
        if selected_names is not None and name not in selected_names:
            continue

        total_count = _population_count(sample_node, f"workspace sample {name!r}")
        sample_subpopulations = _one(
            _direct_children(sample_node, "Subpopulations"),
            f"SampleNode/Subpopulations in workspace sample {name!r}",
        )
        enrichment = _one(
            [
                element
                for element in _direct_children(sample_subpopulations, "Population")
                if element.get("name") == "Human Cell Enrichment"
            ],
            f"direct Human Cell Enrichment population in workspace sample {name!r}",
        )
        enrichment_count = _population_count(
            enrichment, f"Human Cell Enrichment in workspace sample {name!r}"
        )
        if enrichment_count > total_count:
            raise ExtractionError(
                f"Human Cell Enrichment count exceeds total events in workspace sample {name!r}"
            )
        enrichment_subpopulations = _one(
            _direct_children(enrichment, "Subpopulations"),
            f"Human Cell Enrichment/Subpopulations in workspace sample {name!r}",
        )
        human = _one(
            [
                element
                for element in _direct_children(
                    enrichment_subpopulations, "Population"
                )
                if element.get("name") == "HumanCells"
            ],
            f"direct Human Cell Enrichment/HumanCells population in workspace sample {name!r}",
        )
        human_count = _population_count(human, f"HumanCells in workspace sample {name!r}")
        if human_count > enrichment_count:
            raise ExtractionError(
                f"HumanCells count exceeds Human Cell Enrichment count in workspace sample {name!r}"
            )
        subpopulations = _one(
            _direct_children(human, "Subpopulations"),
            f"HumanCells/Subpopulations in workspace sample {name!r}",
        )
        populations: dict[str, int] = {}
        for population in _direct_children(subpopulations, "Population"):
            population_name = (population.get("name") or "").strip()
            if not population_name:
                raise ExtractionError(
                    f"workspace sample {name!r} has an unnamed HumanCells child population"
                )
            if population_name in populations:
                raise ExtractionError(
                    f"workspace sample {name!r} has duplicate direct HumanCells population "
                    f"{population_name!r}"
                )
            populations[population_name] = _population_count(
                population, f"population {population_name!r} in workspace sample {name!r}"
            )
        samples[name] = WorkspaceSample(
            sample_id=sample_id,
            name=name,
            data_set_basename=data_set_basename,
            total_count=total_count,
            human_cell_enrichment_count=enrichment_count,
            human_count=human_count,
            populations=populations,
        )
    if selected_names is not None:
        missing = sorted(selected_names - set(samples))
        if missing:
            raise ExtractionError(
                "workspace lacks crosswalk SampleNode name(s): " + ", ".join(missing)
            )
    return samples


def _percent(numerator: int, denominator: int, context: str) -> float:
    if denominator <= 0:
        raise ExtractionError(f"cannot calculate {context}: denominator is zero")
    return 100.0 * numerator / denominator


def _format_number(value: float) -> str:
    if not math.isfinite(value):
        raise ExtractionError("refusing to write a non-finite result")
    return f"{value:.6f}"


def extract_samples(
    crosswalk: Sequence[CrosswalkRow],
    workspace: Mapping[str, WorkspaceSample],
    fcs_dir: Path,
    expected_cytometer: str,
    expected_dna_channel: str,
    min_human_cells: int,
) -> list[dict[str, object]]:
    results: list[dict[str, object]] = []
    for row in crosswalk:
        if row.wsp_sample_name not in workspace:
            raise ExtractionError(
                f"crosswalk sample {row.mouse_id!r}: workspace has no SampleNode "
                f"named {row.wsp_sample_name!r}"
            )
        wsp = workspace[row.wsp_sample_name]
        fcs_path = fcs_dir / row.fcs_file
        if not fcs_path.is_file():
            raise ExtractionError(
                f"crosswalk sample {row.mouse_id!r}: FCS file does not exist: {fcs_path}"
            )
        fcs = read_fcs_metadata(fcs_path)
        try:
            fcs_acquisition_date = datetime.strptime(
                fcs.keywords["$DATE"], "%d-%b-%Y"
            ).date().isoformat()
        except ValueError as exc:
            raise ExtractionError(
                f"crosswalk sample {row.mouse_id!r}: FCS $DATE "
                f"{fcs.keywords['$DATE']!r} is not DD-MON-YYYY"
            ) from exc
        if fcs_acquisition_date != row.acquisition_date:
            raise ExtractionError(
                f"crosswalk sample {row.mouse_id!r}: acquisition_date "
                f"{row.acquisition_date!r} does not match FCS $DATE "
                f"{fcs_acquisition_date!r}"
            )
        if fcs.keywords["$CYT"] != expected_cytometer:
            raise ExtractionError(
                f"crosswalk sample {row.mouse_id!r}: FCS $CYT {fcs.keywords['$CYT']!r} "
                f"does not equal expected cytometer {expected_cytometer!r}"
            )
        parameter_count = int(fcs.keywords["$PAR"])
        channel_matches = [
            index
            for index in range(1, parameter_count + 1)
            if fcs.keywords.get(f"$P{index}N") == expected_dna_channel
        ]
        if len(channel_matches) != 1:
            raise ExtractionError(
                f"crosswalk sample {row.mouse_id!r}: expected exactly one FCS parameter "
                f"named {expected_dna_channel!r}; found {len(channel_matches)}"
            )
        dna_parameter_index = channel_matches[0]
        fcs_basename = Path(row.fcs_file).name
        if fcs_basename != row.wsp_sample_name:
            raise ExtractionError(
                f"crosswalk sample {row.mouse_id!r}: fcs_file basename {fcs_basename!r} "
                f"does not equal wsp_sample_name {row.wsp_sample_name!r}"
            )
        if wsp.data_set_basename != row.wsp_sample_name:
            raise ExtractionError(
                f"crosswalk sample {row.mouse_id!r}: workspace DataSet basename "
                f"{wsp.data_set_basename!r} does not equal SampleNode name {wsp.name!r}"
            )
        if Path(fcs.keywords["$FIL"]).name != fcs_basename:
            raise ExtractionError(
                f"crosswalk sample {row.mouse_id!r}: FCS $FIL "
                f"{fcs.keywords['$FIL']!r} does not match {fcs_basename!r}"
            )
        fcs_total = int(fcs.keywords["$TOT"])
        if fcs_total != wsp.total_count:
            raise ExtractionError(
                f"crosswalk sample {row.mouse_id!r}: FCS $TOT ({fcs_total}) does not "
                f"equal workspace SampleNode count ({wsp.total_count})"
            )
        required = (*REQUIRED_GATE_NAMES, row.peak_gate_name)
        missing = [name for name in required if name not in wsp.populations]
        if missing:
            raise ExtractionError(
                f"crosswalk sample {row.mouse_id!r}: missing direct HumanCells gate(s): "
                + ", ".join(repr(name) for name in missing)
            )
        for gate_name in set(required):
            if wsp.populations[gate_name] > wsp.human_count:
                raise ExtractionError(
                    f"crosswalk sample {row.mouse_id!r}: gate {gate_name!r} count "
                    "exceeds HumanCells count"
                )

        peak_count = wsp.populations[row.peak_gate_name]
        two_n = wsp.populations["2N"]
        four_n = wsp.populations["4N"]
        mch_two_n = wsp.populations[r"mCh+\2N"]
        mch_four_n = wsp.populations[r"mCh+\4N"]
        summed_mch_gate_counts = mch_two_n + mch_four_n
        human_cells_qc_pass = wsp.human_count >= min_human_cells
        result: dict[str, object] = {
            "mouse_id": row.mouse_id,
            "injected_origin": row.injected_origin,
            "dose_mg_kg": row.dose_mg_kg,
            "fcs_file": row.fcs_file,
            "wsp_sample_name": row.wsp_sample_name,
            "wsp_sample_id": wsp.sample_id,
            "fcs_total_count": fcs_total,
            "wsp_total_count": wsp.total_count,
            "human_cell_enrichment_count": wsp.human_cell_enrichment_count,
            "human_cell_enrichment_pct_total": _format_number(
                _percent(
                    wsp.human_cell_enrichment_count,
                    wsp.total_count,
                    "Human Cell Enrichment percentage",
                )
            ),
            "human_cells_count": wsp.human_count,
            "human_cells_pct_enrichment": _format_number(
                _percent(
                    wsp.human_count,
                    wsp.human_cell_enrichment_count,
                    "HumanCells percentage of Human Cell Enrichment",
                )
            ),
            "human_cells_pct_total": _format_number(
                _percent(wsp.human_count, wsp.total_count, "HumanCells percentage")
            ),
            "human_cells_qc_pass": "TRUE" if human_cells_qc_pass else "FALSE",
            "human_cells_qc_note": (
                "pass"
                if human_cells_qc_pass
                else f"human_cells_count<{min_human_cells}"
            ),
            "peak_gate_name": row.peak_gate_name,
            "reviewed_peak_annotation_n": _format_number(
                float(row.reviewed_peak_annotation_n)
            ),
            "peak_count": peak_count,
            "peak_pct_human_cells": _format_number(
                _percent(peak_count, wsp.human_count, "peak percentage")
            ),
            "two_n_count": two_n,
            "two_n_pct_human_cells": _format_number(
                _percent(two_n, wsp.human_count, "2N percentage")
            ),
            "four_n_count": four_n,
            "four_n_pct_human_cells": _format_number(
                _percent(four_n, wsp.human_count, "4N percentage")
            ),
            "mch_pos_two_n_count": mch_two_n,
            "mch_pos_two_n_pct_human_cells": _format_number(
                _percent(mch_two_n, wsp.human_count, "mCherry-positive 2N percentage")
            ),
            "mch_pos_four_n_count": mch_four_n,
            "mch_pos_four_n_pct_human_cells": _format_number(
                _percent(mch_four_n, wsp.human_count, "mCherry-positive 4N percentage")
            ),
            "summed_mch_pos_two_n_four_n_gate_count": summed_mch_gate_counts,
            "mch_pos_two_n_pct_of_summed_two_n_four_n_gate_counts": _format_number(
                _percent(
                    mch_two_n,
                    summed_mch_gate_counts,
                    "mCherry-positive 2N percentage of summed 2N/4N gate counts",
                )
            ),
            "mch_pos_four_n_pct_of_summed_two_n_four_n_gate_counts": _format_number(
                _percent(
                    mch_four_n,
                    summed_mch_gate_counts,
                    "mCherry-positive 4N percentage of summed 2N/4N gate counts",
                )
            ),
            "fcs_version": fcs.version,
            "fcs_cytometer": fcs.keywords["$CYT"],
            "dna_content_parameter_index": dna_parameter_index,
            "dna_content_parameter_name": expected_dna_channel,
            "crosswalk_acquisition_date": row.acquisition_date,
            "fcs_date": fcs.keywords.get("$DATE", ""),
            "fcs_begin_time": fcs.keywords.get("$BTIM", ""),
            "fcs_end_time": fcs.keywords.get("$ETIM", ""),
        }
        results.append(result)
    return sorted(results, key=lambda item: str(item["mouse_id"]))


def _as_float(row: Mapping[str, object], field: str) -> float:
    return float(str(row[field]))


def _mean_field_or_na(rows: Sequence[Mapping[str, object]], field: str) -> str:
    if not rows:
        return "NA"
    return _format_number(mean(_as_float(row, field) for row in rows))


def _summarize_group(
    rows: Sequence[Mapping[str, object]],
    group_level: str,
    injected_origin: str,
    dose_mg_kg: str,
) -> dict[str, object]:
    qc_rows = [row for row in rows if row["human_cells_qc_pass"] == "TRUE"]
    mch_four_n = sum(int(row["mch_pos_four_n_count"]) for row in rows)
    summed_mch_gate_counts = sum(
        int(row["summed_mch_pos_two_n_four_n_gate_count"]) for row in rows
    )
    return {
        "group_level": group_level,
        "injected_origin": injected_origin,
        "dose_mg_kg": dose_mg_kg,
        "n_samples": len(rows),
        "n_samples_human_cells_qc_pass": len(qc_rows),
        "total_human_cells_count": sum(int(row["human_cells_count"]) for row in rows),
        "total_summed_mch_pos_two_n_four_n_gate_count": summed_mch_gate_counts,
        "reviewed_peak_gate_labels": ";".join(
            sorted({str(row["peak_gate_name"]) for row in rows})
        ),
        "mean_peak_pct_human_cells_all_samples": _format_number(
            mean(_as_float(row, "peak_pct_human_cells") for row in rows)
        ),
        "mean_peak_pct_human_cells_qc_pass": _mean_field_or_na(
            qc_rows, "peak_pct_human_cells"
        ),
        "mean_two_n_pct_human_cells": _format_number(
            mean(_as_float(row, "two_n_pct_human_cells") for row in rows)
        ),
        "mean_four_n_pct_human_cells": _format_number(
            mean(_as_float(row, "four_n_pct_human_cells") for row in rows)
        ),
        "mean_mch_pos_four_n_pct_of_summed_two_n_four_n_gate_counts": _format_number(
            mean(
                _as_float(
                    row,
                    "mch_pos_four_n_pct_of_summed_two_n_four_n_gate_counts",
                )
                for row in rows
            )
        ),
        "pooled_mch_pos_four_n_pct_of_summed_two_n_four_n_gate_counts": _format_number(
            _percent(
                mch_four_n,
                summed_mch_gate_counts,
                "pooled mCherry-positive 4N percentage of summed 2N/4N gate counts",
            )
        ),
    }


def summarize_groups(rows: Sequence[Mapping[str, object]]) -> list[dict[str, object]]:
    if not rows:
        raise ExtractionError("cannot summarize an empty sample table")
    summaries = [_summarize_group(rows, "all", "all", "all")]
    origins = sorted({str(row["injected_origin"]) for row in rows})
    for origin in origins:
        origin_rows = [row for row in rows if row["injected_origin"] == origin]
        summaries.append(_summarize_group(origin_rows, "injected_origin", origin, "all"))
        doses = sorted(
            {str(row["dose_mg_kg"]) for row in origin_rows}, key=lambda value: float(value)
        )
        for dose in doses:
            dose_rows = [row for row in origin_rows if str(row["dose_mg_kg"]) == dose]
            summaries.append(
                _summarize_group(dose_rows, "injected_origin_x_dose", origin, dose)
            )
    return summaries


def extract_paired_sensitivity(
    sensitivity: Sequence[SensitivityRow],
    crosswalk: Sequence[CrosswalkRow],
    workspace: Mapping[str, WorkspaceSample],
) -> list[dict[str, object]]:
    crosswalk_by_mouse = {row.mouse_id: row for row in crosswalk}
    results = []
    for row in sensitivity:
        selected_contract = crosswalk_by_mouse.get(row.mouse_id)
        if selected_contract is None:
            raise ExtractionError(
                f"paired-sensitivity mouse {row.mouse_id!r} is absent from the crosswalk"
            )
        if selected_contract.injected_origin != "4N":
            raise ExtractionError(
                f"paired-sensitivity mouse {row.mouse_id!r} is not 4N-origin"
            )
        selected_expected = (
            selected_contract.wsp_sample_name,
            selected_contract.peak_gate_name,
            float(selected_contract.reviewed_peak_annotation_n),
        )
        selected_observed = (
            row.selected_wsp_sample_name,
            row.selected_peak_gate_name,
            float(row.selected_peak_annotation_n),
        )
        if selected_observed != selected_expected:
            raise ExtractionError(
                f"paired-sensitivity selected mapping disagrees with crosswalk for "
                f"{row.mouse_id!r}"
            )
        if row.paired_wsp_sample_name == row.selected_wsp_sample_name:
            raise ExtractionError(
                f"paired-sensitivity sample is not distinct for {row.mouse_id!r}"
            )
        selected = workspace.get(row.selected_wsp_sample_name)
        paired = workspace.get(row.paired_wsp_sample_name)
        if selected is None or paired is None:
            raise ExtractionError(
                f"paired-sensitivity workspace mapping is incomplete for {row.mouse_id!r}"
            )
        for sample, gate_name, label in (
            (selected, row.selected_peak_gate_name, "selected"),
            (paired, row.paired_peak_gate_name, "paired"),
        ):
            if gate_name not in sample.populations:
                raise ExtractionError(
                    f"paired-sensitivity {label} gate {gate_name!r} is absent for "
                    f"{row.mouse_id!r}"
                )
            if sample.populations[gate_name] > sample.human_count:
                raise ExtractionError(
                    f"paired-sensitivity {label} gate count exceeds HumanCells for "
                    f"{row.mouse_id!r}"
                )
        selected_peak_count = selected.populations[row.selected_peak_gate_name]
        paired_peak_count = paired.populations[row.paired_peak_gate_name]
        results.append(
            {
                "mouse_id": row.mouse_id,
                "selected_wsp_sample_name": row.selected_wsp_sample_name,
                "selected_peak_gate_name": row.selected_peak_gate_name,
                "selected_peak_annotation_n": _format_number(
                    float(row.selected_peak_annotation_n)
                ),
                "selected_human_cells_count": selected.human_count,
                "selected_peak_count": selected_peak_count,
                "selected_peak_pct_human_cells": _format_number(
                    _percent(
                        selected_peak_count,
                        selected.human_count,
                        "selected sensitivity peak percentage",
                    )
                ),
                "paired_wsp_sample_name": row.paired_wsp_sample_name,
                "paired_peak_gate_name": row.paired_peak_gate_name,
                "paired_peak_annotation_n": _format_number(
                    float(row.paired_peak_annotation_n)
                ),
                "paired_human_cells_count": paired.human_count,
                "paired_peak_count": paired_peak_count,
                "paired_peak_pct_human_cells": _format_number(
                    _percent(
                        paired_peak_count,
                        paired.human_count,
                        "paired sensitivity peak percentage",
                    )
                ),
                "relationship_status": row.relationship_status,
                "paired_fcs_imported": "FALSE",
            }
        )
    return sorted(results, key=lambda item: str(item["mouse_id"]))


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    try:
        with path.open("rb") as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(chunk)
    except OSError as exc:
        raise ExtractionError(f"cannot hash input {path}: {exc}") from exc
    return digest.hexdigest()


def verify_expected_input_manifest(
    manifest_path: Path,
    crosswalk_path: Path,
    sensitivity_path: Path,
    workspace_path: Path,
    crosswalk: Sequence[CrosswalkRow],
    fcs_dir: Path,
) -> list[dict[str, object]]:
    try:
        with manifest_path.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if tuple(reader.fieldnames or ()) != EXPECTED_MANIFEST_COLUMNS:
                raise ExtractionError(
                    "expected-input manifest columns must be exactly: "
                    + ", ".join(EXPECTED_MANIFEST_COLUMNS)
                )
            manifest_rows = list(reader)
    except OSError as exc:
        raise ExtractionError(
            f"cannot read expected-input manifest {manifest_path}: {exc}"
        ) from exc

    manifest_root = manifest_path.parent.resolve()
    expected: dict[tuple[str, Path], dict[str, object]] = {}
    for line_number, raw in enumerate(manifest_rows, start=2):
        role = (raw.get("input_role") or "").strip()
        relative_text = (raw.get("relative_path") or "").strip()
        size_text = (raw.get("size_bytes") or "").strip()
        digest = (raw.get("sha256") or "").strip().lower()
        if role not in {"crosswalk", "selection_sensitivity", "workspace", "fcs"}:
            raise ExtractionError(
                f"expected-input manifest line {line_number}: invalid input_role {role!r}"
            )
        relative_path = Path(relative_text)
        if (
            not relative_text
            or relative_path.is_absolute()
            or ".." in relative_path.parts
        ):
            raise ExtractionError(
                f"expected-input manifest line {line_number}: relative_path must stay "
                "within the dataset directory"
            )
        resolved = (manifest_root / relative_path).resolve()
        try:
            resolved.relative_to(manifest_root)
        except ValueError as exc:
            raise ExtractionError(
                f"expected-input manifest line {line_number}: resolved path escapes "
                "the dataset directory"
            ) from exc
        try:
            size = int(size_text)
        except ValueError as exc:
            raise ExtractionError(
                f"expected-input manifest line {line_number}: size_bytes is not an integer"
            ) from exc
        if size < 0 or not re.fullmatch(r"[0-9a-f]{64}", digest):
            raise ExtractionError(
                f"expected-input manifest line {line_number}: invalid size or SHA-256"
            )
        key = (role, resolved)
        if key in expected:
            raise ExtractionError(
                f"expected-input manifest has duplicate role/path at line {line_number}"
            )
        expected[key] = {
            "input_role": role,
            "logical_path": relative_path.as_posix(),
            "size_bytes": size,
            "sha256": digest,
        }

    observed_paths: list[tuple[str, Path]] = [
        ("crosswalk", crosswalk_path.resolve()),
        ("selection_sensitivity", sensitivity_path.resolve()),
        ("workspace", workspace_path.resolve()),
    ]
    observed_paths.extend(
        ("fcs", (fcs_dir / row.fcs_file).resolve()) for row in crosswalk
    )
    if len(expected) != len(observed_paths):
        raise ExtractionError(
            f"expected-input manifest has {len(expected)} rows; expected exactly "
            f"{len(observed_paths)} for the selected cohort"
        )
    missing = sorted(
        f"{role}:{path.name}"
        for role, path in observed_paths
        if (role, path) not in expected
    )
    unexpected = sorted(
        f"{role}:{path.name}"
        for role, path in expected
        if (role, path) not in set(observed_paths)
    )
    if missing or unexpected:
        raise ExtractionError(
            "expected-input manifest does not exactly match selected inputs: "
            f"missing={missing}; unexpected={unexpected}"
        )

    rows = []
    for role, path in observed_paths:
        recorded = expected[(role, path)]
        try:
            observed_size = path.stat().st_size
        except OSError as exc:
            raise ExtractionError(f"cannot stat input {path}: {exc}") from exc
        observed_digest = sha256_file(path)
        if (
            observed_size != recorded["size_bytes"]
            or observed_digest != recorded["sha256"]
        ):
            raise ExtractionError(
                f"input does not match expected size/SHA-256 manifest: "
                f"{recorded['logical_path']}"
            )
        rows.append(
            {
                "input_role": role,
                "logical_path": recorded["logical_path"],
                "size_bytes": observed_size,
                "sha256": observed_digest,
            }
        )
    role_order = {
        "crosswalk": 0,
        "selection_sensitivity": 1,
        "workspace": 2,
        "fcs": 3,
    }
    return sorted(
        rows,
        key=lambda row: (
            role_order[str(row["input_role"])],
            str(row["logical_path"]),
        ),
    )


def _write_tsv_atomic(
    path: Path, fieldnames: Sequence[str], rows: Iterable[Mapping[str, object]]
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=path.parent
    )
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=fieldnames,
                delimiter="\t",
                lineterminator="\n",
                extrasaction="raise",
            )
            writer.writeheader()
            writer.writerows(rows)
        os.replace(temporary_name, path)
    except Exception:
        try:
            os.unlink(temporary_name)
        except FileNotFoundError:
            pass
        raise


def run(args: argparse.Namespace) -> tuple[Path, Path, Path, Path]:
    crosswalk_path = Path(args.crosswalk).resolve()
    sensitivity_path = Path(args.paired_sensitivity).resolve()
    workspace_path = Path(args.workspace).resolve()
    fcs_dir = Path(args.fcs_dir).resolve()
    expected_manifest_path = Path(args.expected_input_manifest).resolve()
    output_dir = Path(args.output_dir).resolve()
    if args.expected_samples <= 0:
        raise ExtractionError("--expected-samples must be a positive integer")
    if not workspace_path.is_file():
        raise ExtractionError(f"workspace does not exist: {workspace_path}")
    if not sensitivity_path.is_file():
        raise ExtractionError(
            f"paired-sensitivity table does not exist: {sensitivity_path}"
        )
    if not fcs_dir.is_dir():
        raise ExtractionError(f"FCS directory does not exist: {fcs_dir}")
    if not expected_manifest_path.is_file():
        raise ExtractionError(
            f"expected-input manifest does not exist: {expected_manifest_path}"
        )

    crosswalk = read_crosswalk(crosswalk_path, args.expected_samples)
    sensitivity = read_sensitivity_table(sensitivity_path)
    hashes = verify_expected_input_manifest(
        expected_manifest_path,
        crosswalk_path,
        sensitivity_path,
        workspace_path,
        crosswalk,
        fcs_dir,
    )
    workspace = read_workspace(
        workspace_path,
        {row.wsp_sample_name for row in crosswalk}
        | {row.paired_wsp_sample_name for row in sensitivity},
    )
    if not args.expected_cytometer.strip():
        raise ExtractionError("--expected-cytometer cannot be empty")
    if not args.expected_dna_channel.strip():
        raise ExtractionError("--expected-dna-channel cannot be empty")
    if args.min_human_cells < 0:
        raise ExtractionError("--min-human-cells must be nonnegative")
    samples = extract_samples(
        crosswalk,
        workspace,
        fcs_dir,
        expected_cytometer=args.expected_cytometer,
        expected_dna_channel=args.expected_dna_channel,
        min_human_cells=args.min_human_cells,
    )
    summaries = summarize_groups(samples)
    sensitivity_results = extract_paired_sensitivity(
        sensitivity, crosswalk, workspace
    )

    sample_path = output_dir / "endpoint_flow_per_sample.tsv"
    summary_path = output_dir / "endpoint_flow_group_summary.tsv"
    hashes_path = output_dir / "input_hashes.tsv"
    sensitivity_output_path = output_dir / "paired_workspace_sensitivity.tsv"
    _write_tsv_atomic(sample_path, PER_SAMPLE_FIELDS, samples)
    _write_tsv_atomic(summary_path, GROUP_FIELDS, summaries)
    _write_tsv_atomic(
        hashes_path,
        ("input_role", "logical_path", "size_bytes", "sha256"),
        hashes,
    )
    _write_tsv_atomic(
        sensitivity_output_path,
        SENSITIVITY_OUTPUT_FIELDS,
        sensitivity_results,
    )
    return sample_path, summary_path, hashes_path, sensitivity_output_path


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Extract reviewed endpoint DNA-content gate counts from an explicit "
            "16-sample FlowJo/FCS crosswalk."
        )
    )
    parser.add_argument("--crosswalk", required=True, help="Explicit sample crosswalk TSV")
    parser.add_argument(
        "--paired-sensitivity",
        required=True,
        help="Explicit unresolved paired-acquisition sensitivity mapping",
    )
    parser.add_argument("--workspace", required=True, help="Reviewed FlowJo .wsp file")
    parser.add_argument(
        "--expected-input-manifest",
        required=True,
        help="Pinned dataset-relative size/SHA-256 manifest for all selected inputs",
    )
    parser.add_argument(
        "--fcs-dir", required=True, help="Directory containing crosswalk-relative FCS files"
    )
    parser.add_argument("--output-dir", required=True, help="Directory for deterministic TSVs")
    parser.add_argument(
        "--expected-samples",
        type=int,
        default=16,
        help="Required crosswalk row count (default: 16)",
    )
    parser.add_argument(
        "--expected-cytometer",
        default="FACSCantoII",
        help="Required FCS $CYT value (default: FACSCantoII)",
    )
    parser.add_argument(
        "--expected-dna-channel",
        default="450/50 Violet B-A",
        help="Required unique FCS parameter name for DNA content",
    )
    parser.add_argument(
        "--min-human-cells",
        type=int,
        default=1000,
        help="Flag, but do not exclude, samples below this HumanCells count (default: 1000)",
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        paths = run(args)
    except ExtractionError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    for path in paths:
        print(f"Wrote {path.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
