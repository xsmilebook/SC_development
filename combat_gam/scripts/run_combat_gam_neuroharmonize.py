#!/usr/bin/env python3
import argparse
import os
import sys
import numpy as np
import pandas as pd
import pyreadr

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
    data = pyreadr.read_r(path)
    if not data:
        raise ValueError(f"No objects found in {path}")
    return next(iter(data.values()))


def main():
    parser = argparse.ArgumentParser(description="Run ComBat-GAM using neuroHarmonize on SC tables.")
    parser.add_argument("--input-rds", required=True)
    parser.add_argument("--output-rds", required=True)
    parser.add_argument("--batch-col", required=True)
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

    covar_cols = [args.id_col, args.batch_col, args.age_col, args.sex_col, args.mean_fd_col]
    extra_cols = [c for c in args.extra_covars.split(",") if c]
    covar_cols.extend(extra_cols)
    covar_cols = list(dict.fromkeys(covar_cols))

    needed_cols = sc_cols + covar_cols
    missing = [c for c in needed_cols if c not in df.columns]
    if missing:
        raise ValueError(f"Missing columns in input: {missing}")

    df = df[needed_cols].dropna()
    if args.test_n > 0:
        df = df.sample(n=min(args.test_n, len(df)), random_state=args.seed)

    covars = pd.DataFrame()
    covars["SITE"] = df[args.batch_col].astype(str)
    covars[args.age_col] = _ensure_numeric(df[args.age_col])
    covars[args.sex_col] = _ensure_numeric(df[args.sex_col])
    covars[args.mean_fd_col] = _ensure_numeric(df[args.mean_fd_col])
    for col in extra_cols:
        covars[col] = _ensure_numeric(df[col])

    data_matrix = df[sc_cols].to_numpy()

    _, data_adj = harmonizationLearn(
        data_matrix,
        covars,
        smooth_terms=[args.age_col],
        smooth_df=4,
        smooth_degree=3,
        smooth_fx=True,
    )

    out = df[[args.id_col]].copy()
    for col in [args.batch_col, args.age_col, args.sex_col, args.mean_fd_col] + extra_cols:
        out[col] = df[col].values
    for idx, col in enumerate(sc_cols):
        out[f"{col}_h"] = data_adj[:, idx]

    pyreadr.write_rds(args.output_rds, out)


if __name__ == "__main__":
    main()
