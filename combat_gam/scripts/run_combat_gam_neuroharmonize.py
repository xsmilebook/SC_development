#!/usr/bin/env python3
import argparse
import os
import sys
import numpy as np
import pandas as pd
import pyreadr
import rpy2.robjects as ro
from rpy2.robjects import default_converter, pandas2ri
from rpy2.robjects.conversion import localconverter
from rpy2.robjects.packages import importr

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, "..", ".."))
NH_ROOT = os.path.join(REPO_ROOT, "combat_gam", "neuroHarmonize")
if NH_ROOT not in sys.path:
    sys.path.insert(0, NH_ROOT)

from neuroHarmonize import harmonizationLearn


def _ensure_numeric(series):
    if pd.api.types.is_numeric_dtype(series):
        return series.astype(float)
    return pd.Categorical(series).codes.astype(float)


def _load_rds(path):
    """
    Load an RDS file into a pandas.DataFrame.

    Primary backend: pyreadr (fast).
    Fallback: R's readRDS via rpy2 (handles newer RDS features that librdata/pyreadr may not support).
    """
    try:
        data = pyreadr.read_r(path)
        if not data:
            raise ValueError(f"No objects found in {path}")
        obj = next(iter(data.values()))
        if not isinstance(obj, pd.DataFrame):
            raise ValueError(f"Expected data.frame in {path}, got {type(obj)}")
        return obj
    except Exception as e:
        # Common failure on HPC: pyreadr/librdata can't parse RDS created with newer R versions
        # (e.g., unsupported serialization features). Fall back to R itself.
        sys.stderr.write(f"[WARN] pyreadr failed to read {path}: {e}\n")
        sys.stderr.write("[WARN] Falling back to rpy2 + readRDS(). This is slower but more compatible.\n")
        try:
            ro.globalenv[".rds_path"] = ro.StrVector([path])[0]
            r_obj = ro.r('as.data.frame(readRDS(.rds_path))')
            with localconverter(default_converter + pandas2ri.converter):
                df = ro.conversion.rpy2py(r_obj)
            if not isinstance(df, pd.DataFrame):
                raise ValueError(f"rpy2 readRDS did not return a DataFrame for {path}: {type(df)}")
            return df
        except Exception as e2:
            raise RuntimeError(f"Failed to read RDS via both pyreadr and rpy2: {path}. pyreadr_error={e} rpy2_error={e2}") from e2


def _build_mgcv_basis(age_values, age_col):
    importr("mgcv")
    age_values = np.asarray(age_values, dtype=float)
    with localconverter(default_converter + pandas2ri.converter):
        r_df = ro.conversion.py2rpy(pd.DataFrame({age_col: age_values}))
    ro.globalenv["df"] = r_df
    r_expr = f'sm <- smoothCon(s(`{age_col}`, k=3, bs="tp", fx=TRUE), data=df)[[1]]; sm$X'
    x_mat = np.asarray(ro.r(r_expr))
    stds = x_mat.std(axis=0)
    keep = stds > 0
    x_mat = x_mat[:, keep]
    columns = [f"{age_col}_s_tp_{idx + 1}" for idx in range(x_mat.shape[1])]
    return pd.DataFrame(x_mat, columns=columns)


def _sample_test_df(df, batch_col, sex_col, test_n, seed):
    if test_n <= 0 or test_n >= len(df):
        return df
    batches = df[batch_col].dropna().unique()
    if len(batches) < 2:
        return df.sample(n=min(test_n, len(df)), random_state=seed)
    n_per = max(2, test_n // len(batches))
    pieces = []
    for batch in batches:
        sub = df[df[batch_col] == batch]
        n_take = min(n_per, len(sub))
        pieces.append(sub.sample(n=n_take, random_state=seed))
    sampled = pd.concat(pieces)
    remaining = df.drop(sampled.index)
    if len(sampled) < test_n and not remaining.empty:
        extra = remaining.sample(n=min(test_n - len(sampled), len(remaining)), random_state=seed)
        sampled = pd.concat([sampled, extra])
    if sampled[sex_col].nunique() < 2 or sampled[batch_col].nunique() < 2:
        sampled = df.sample(n=min(test_n, len(df)), random_state=seed)
    if sampled[sex_col].nunique() < 2 or sampled[batch_col].nunique() < 2:
        return df
    return sampled


def main():
    parser = argparse.ArgumentParser(description="Run ComBat-GAM using neuroHarmonize on SC tables.")
    parser.add_argument("--input-rds", required=True)
    parser.add_argument("--output-rds", required=True)
    parser.add_argument("--batch-col", required=True)
    parser.add_argument(
        "--batch-source-rds",
        default="",
        help="Optional RDS providing the batch column (e.g., site). If input lacks --batch-col, "
        "we will merge it from this table using --batch-source-id-col and --batch-source-batch-col.",
    )
    parser.add_argument("--batch-source-id-col", default="", help="ID column in --batch-source-rds (defaults to --id-col).")
    parser.add_argument(
        "--batch-source-batch-col", default="", help="Batch column in --batch-source-rds (defaults to --batch-col)."
    )
    parser.add_argument("--id-col", required=True)
    parser.add_argument("--age-col", required=True)
    parser.add_argument("--sex-col", required=True)
    parser.add_argument("--mean-fd-col", required=True)
    parser.add_argument("--extra-covars", default="")
    parser.add_argument("--test-n", type=int, default=0)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    df = _load_rds(args.input_rds)
    sc_cols = [c for c in df.columns if c.startswith("SC.")]
    if not sc_cols:
        raise ValueError("No SC.* columns found in input data.")

    # Optional: inject missing batch column from a "source" RDS (common for derived tables).
    if args.batch_col not in df.columns and args.batch_source_rds:
        src_id_col = args.batch_source_id_col or args.id_col
        src_batch_col = args.batch_source_batch_col or args.batch_col
        src = _load_rds(args.batch_source_rds)
        if src_id_col not in src.columns or src_batch_col not in src.columns:
            raise ValueError(
                f"batch-source-rds missing required columns: {src_id_col}, {src_batch_col}. "
                f"Available: {list(src.columns)}"
            )
        src_map = src[[src_id_col, src_batch_col]].dropna().drop_duplicates()
        before = df.shape[0]
        df = df.merge(src_map, how="left", left_on=args.id_col, right_on=src_id_col)
        # If id column name differs, drop the extra join key.
        if src_id_col != args.id_col and src_id_col in df.columns:
            df = df.drop(columns=[src_id_col])
        if args.batch_col not in df.columns:
            # If the source batch column name is different, rename it to the requested batch-col.
            if src_batch_col in df.columns:
                df = df.rename(columns={src_batch_col: args.batch_col})
        missing_n = int(df[args.batch_col].isna().sum()) if args.batch_col in df.columns else before
        sys.stderr.write(
            f"[INFO] Injected batch column '{args.batch_col}' from {args.batch_source_rds}; "
            f"rows={before}, missing_batch_after_merge={missing_n}\n"
        )
        if args.batch_col not in df.columns or missing_n == before:
            raise ValueError(
                f"Failed to inject batch column '{args.batch_col}' from {args.batch_source_rds}. "
                f"Check ID alignment between {args.id_col} (input) and {src_id_col} (source)."
            )

    covar_cols = [args.id_col, args.batch_col, args.age_col, args.sex_col, args.mean_fd_col]
    extra_cols = [c for c in args.extra_covars.split(",") if c]
    covar_cols.extend(extra_cols)
    covar_cols = list(dict.fromkeys(covar_cols))

    needed_cols = sc_cols + covar_cols
    missing = [c for c in needed_cols if c not in df.columns]
    if missing:
        raise ValueError(f"Missing columns in input: {missing}")

    df = df[needed_cols].dropna()
    df = _sample_test_df(df, args.batch_col, args.sex_col, args.test_n, args.seed)

    covars = pd.DataFrame()
    covars["SITE"] = df[args.batch_col].astype(str)
    covars[args.sex_col] = _ensure_numeric(df[args.sex_col])
    covars[args.mean_fd_col] = _ensure_numeric(df[args.mean_fd_col])
    for col in extra_cols:
        covars[col] = _ensure_numeric(df[col])

    basis = _build_mgcv_basis(df[args.age_col], args.age_col)
    covars = pd.concat([covars.reset_index(drop=True), basis.reset_index(drop=True)], axis=1)

    data_matrix = df[sc_cols].to_numpy()

    _, data_adj = harmonizationLearn(
        data_matrix,
        covars,
        smooth_terms=[],
    )

    out = df[[args.id_col]].copy()
    for col in [args.batch_col, args.age_col, args.sex_col, args.mean_fd_col] + extra_cols:
        out[col] = df[col].values
    for idx, col in enumerate(sc_cols):
        out[f"{col}_h"] = data_adj[:, idx]

    pyreadr.write_rds(args.output_rds, out)


if __name__ == "__main__":
    main()
