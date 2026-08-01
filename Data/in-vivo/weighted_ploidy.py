#!/usr/bin/env python3
"""
Compute chromosome-length–weighted ploidy per cell from one or more *.sps.cbs files.

Supports two common formats:

(A) "Wide" format (like the provided example):
    - rows: cell IDs
    - columns: segments named like "1:0-1.25e+08" or "7:17298652-152759749"
    - values: integer/float copy number per segment (may contain NA)

(B) "Long" format (used inside calibrate_ploidy.py):
    - columns include: chr, seglength, and one or more cell columns starting with "SP_"
    - each row corresponds to a chromosome (or aggregated chromosome), seglength is chromosome length

Ploidy is computed as:
    sum_i( CN_i * length_i ) / sum_i( length_i )
where i indexes segments/chromosomes and the sums are taken over selected chromosomes (default 1–22).
Missing CN values are ignored on a per-cell basis (denominator only includes lengths for finite CN).

An additional column, total_chromosomes, reports the estimated total number of selected chromosomes
per cell. This is computed by taking a simple (unweighted) mean copy number across segments within
each selected chromosome and then summing those chromosome-level means across chromosomes.

Example:
    python ../weighted_ploidy.py SUM159-2N-30-0_harvest.sps.cbs -o ploidy.csv --hist ploidy_hist.png
"""
from __future__ import annotations

import argparse
import os
import re
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

# Use non-interactive backend for servers / headless environments.
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


_SEG_RE = re.compile(r'^(?P<chr>[^:]+):(?P<start>[^-]+)-(?P<end>.+)$')


def _coerce_chr(ch: str) -> Optional[int]:
    ch = str(ch).strip()
    if ch.lower() == "x":
        return 23
    if ch.lower() == "y":
        return 24
    try:
        return int(ch)
    except Exception:
        return None


def _parse_segment_col(col: str) -> Optional[Tuple[int, float]]:
    """
    Parse a segment column name like '1:0-1.25e+08' into (chr_int, length).
    Returns None if parsing fails.
    """
    m = _SEG_RE.match(str(col))
    if not m:
        return None
    chr_i = _coerce_chr(m.group("chr"))
    if chr_i is None:
        return None
    try:
        start = float(m.group("start"))
        end = float(m.group("end"))
    except Exception:
        return None
    length = float(end - start)
    if not np.isfinite(length) or length <= 0:
        return None
    return chr_i, length


def _select_chr(chr_i: int, include_sex: bool) -> bool:
    if 1 <= chr_i <= 22:
        return True
    if include_sex and chr_i in (23, 24):
        return True
    return False


def _unweighted_chr_totals_from_matrix(
    X: np.ndarray,
    chr_ids: np.ndarray,
    selected_chr_order: List[int],
) -> np.ndarray:
    """
    Estimate total chromosome count per cell by:
      1) computing the simple (unweighted) mean CN within each selected chromosome
      2) summing those chromosome-level means across chromosomes

    Parameters
    ----------
    X : ndarray, shape (n_cells, n_segments)
        Copy-number matrix.
    chr_ids : ndarray, shape (n_segments,)
        Integer chromosome ID for each segment/row of X.
    selected_chr_order : list[int]
        Chromosomes to include in the sum.

    Returns
    -------
    ndarray, shape (n_cells,)
        Estimated total chromosome count per cell.
    """
    totals = np.zeros(X.shape[0], dtype=float)
    for chr_i in selected_chr_order:
        mask = chr_ids == chr_i
        if not np.any(mask):
            continue
        X_chr = X[:, mask]
        mean_chr = np.nanmean(X_chr, axis=1)
        totals += np.nan_to_num(mean_chr, nan=0.0)
    return totals


def ploidy_from_wide_sps_cbs(
    path: str,
    include_sex: bool = False,
    min_frac_covered: float = 0.0,
) -> pd.DataFrame:
    """
    Wide format:
        index: cell IDs
        columns: segment names (chr:start-end)
        values: CN (may contain NA)
    Returns DataFrame with columns: cell_id, ploidy, total_chromosomes, frac_covered
    """
    df = pd.read_csv(path, sep="\t", index_col=0)
    if df.shape[1] == 0:
        raise ValueError(f"No segment columns found in {path}")

    seg_info: Dict[str, Tuple[int, float]] = {}
    for col in df.columns:
        parsed = _parse_segment_col(col)
        if parsed is not None:
            seg_info[col] = parsed

    if not seg_info:
        raise ValueError(f"Could not parse any segment columns like 'chr:start-end' in {path}")

    sel_cols: List[str] = []
    weights: List[float] = []
    chr_ids: List[int] = []
    selected_chr_order: List[int] = []
    seen_chr = set()
    for col, (chr_i, length) in seg_info.items():
        if _select_chr(chr_i, include_sex):
            sel_cols.append(col)
            weights.append(length)
            chr_ids.append(chr_i)
            if chr_i not in seen_chr:
                selected_chr_order.append(chr_i)
                seen_chr.add(chr_i)

    if not sel_cols:
        raise ValueError(f"No usable autosomal (or selected) segments found in {path}")

    w = np.asarray(weights, dtype=float)
    chr_ids_arr = np.asarray(chr_ids, dtype=int)
    X = df[sel_cols].apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float)

    finite = np.isfinite(X)
    denom = (finite * w).sum(axis=1)
    denom_total = float(w.sum())
    frac_covered = denom / denom_total

    numer = np.nansum(X * w, axis=1)
    ploidy = np.full(X.shape[0], np.nan, dtype=float)
    ok = denom > 0
    ploidy[ok] = numer[ok] / denom[ok]

    total_chromosomes = _unweighted_chr_totals_from_matrix(
        X=X,
        chr_ids=chr_ids_arr,
        selected_chr_order=selected_chr_order,
    )

    if min_frac_covered > 0:
        bad = frac_covered < min_frac_covered
        ploidy[bad] = np.nan
        total_chromosomes[bad] = np.nan

    out = pd.DataFrame(
        {
            "cell_id": df.index.astype(str),
            "ploidy": ploidy,
            "total_chromosomes": total_chromosomes,
            "frac_covered": frac_covered,
        }
    )
    return out


def ploidy_from_long_sps_cbs(path: str, include_sex: bool = False) -> pd.DataFrame:
    """
    Long format:
        columns: chr, seglength, SP_* (one per cell)
    Returns DataFrame with columns: cell_id, ploidy, total_chromosomes, frac_covered
    """
    cbs = pd.read_csv(path, sep="\t")
    required = {"chr", "seglength"}
    if not required.issubset(set(cbs.columns)):
        raise ValueError("Not a long-format file (missing chr/seglength columns).")

    cell_cols = [c for c in cbs.columns if str(c).startswith("SP_")]
    if not cell_cols:
        raise ValueError("Not a long-format file (no SP_* columns found).")

    df = cbs[["chr", "seglength"] + cell_cols].copy()
    df["chr"] = df["chr"].apply(_coerce_chr)
    df = df[df["chr"].notna()]
    df["chr"] = df["chr"].astype(int)
    df = df[df["chr"].apply(lambda x: _select_chr(int(x), include_sex))]

    if df.empty:
        raise ValueError("No usable chromosomes remain after filtering (1–22, optionally 23/24).")

    selected_chr_order = sorted(df["chr"].unique().tolist())
    seglen = df["seglength"].astype(float).to_numpy()
    chr_ids = df["chr"].astype(int).to_numpy()
    X = df[cell_cols].apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float)

    numer = np.nansum(X * seglen[:, None], axis=0)
    denom = np.nansum(np.isfinite(X) * seglen[:, None], axis=0)
    Ltot = float(seglen.sum())

    ploidy = np.full(X.shape[1], np.nan, dtype=float)
    ok = denom > 0
    ploidy[ok] = numer[ok] / denom[ok]

    total_chromosomes = _unweighted_chr_totals_from_matrix(
        X=X.T,
        chr_ids=chr_ids,
        selected_chr_order=selected_chr_order,
    )

    out = pd.DataFrame(
        {
            "cell_id": cell_cols,
            "ploidy": ploidy,
            "total_chromosomes": total_chromosomes,
            "frac_covered": denom / Ltot,
        }
    )
    return out


def ploidy_from_file(
    path: str,
    include_sex: bool = False,
    min_frac_covered: float = 0.0,
) -> pd.DataFrame:
    """
    Detect file format and compute ploidy per cell.
    """
    try:
        out = ploidy_from_wide_sps_cbs(path, include_sex=include_sex, min_frac_covered=min_frac_covered)
        out["format"] = "wide"
        return out
    except Exception:
        out = ploidy_from_long_sps_cbs(path, include_sex=include_sex)
        out["format"] = "long"
        return out


def plot_histogram(values: np.ndarray, out_png: str, bins: int = 60, title: str = "Ploidy distribution"):
    v = np.asarray(values, dtype=float)
    v = v[np.isfinite(v)]
    plt.figure(figsize=(6, 4))
    plt.hist(v, bins=bins)
    plt.xlabel("chr-length weighted ploidy")
    plt.ylabel("count")
    plt.title(title)
    plt.tight_layout()
    os.makedirs(os.path.dirname(out_png) or ".", exist_ok=True)
    plt.savefig(out_png, dpi=200)
    plt.close()


def main(argv: Optional[List[str]] = None) -> int:
    ap = argparse.ArgumentParser(description="Compute chromosome-length–weighted ploidy per cell from *.sps.cbs files.")
    ap.add_argument("files", nargs="+", help="One or more input files (e.g., *.sps.cbs).")
    ap.add_argument("-o", "--out", default=None, help="Output CSV/TSV path (default: print TSV to stdout).")
    ap.add_argument("--sep", default="\t", choices=["\t", ",", "tsv", "csv"], help="Output separator.")
    ap.add_argument("--include-sex", action="store_true", help="Include chr X/Y (23/24) in weighting.")
    ap.add_argument(
        "--min-frac-covered",
        type=float,
        default=0.0,
        help="Set ploidy and total_chromosomes to NA if the covered fraction of genome length is below this threshold (0–1).",
    )
    ap.add_argument("--hist", default=None, help="If set, save a histogram PNG of ploidy (all files combined).")
    ap.add_argument("--bins", type=int, default=60, help="Histogram bins (if --hist is provided).")
    args = ap.parse_args(argv)

    sep = args.sep
    if sep == "tsv":
        sep = "\t"
    elif sep == "csv":
        sep = ","

    rows: List[pd.DataFrame] = []
    for f in args.files:
        df = ploidy_from_file(f, include_sex=args.include_sex, min_frac_covered=args.min_frac_covered)
        df.insert(0, "file", os.path.basename(f))
        rows.append(df)

    out = pd.concat(rows, ignore_index=True)

    if args.hist:
        plot_histogram(out["ploidy"].to_numpy(), args.hist, bins=args.bins, title="Ploidy distribution (chr-length weighted)")

    if args.out:
        out.to_csv(args.out, sep=sep, index=False)
    else:
        print(out.to_csv(sep="\t", index=False), end="")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
