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

from neuroHarmonize import harmonizationLearn  # noqa: E402


def _ensure_numeric(series: pd.Series) -> np.ndarray:
    if pd.api.types.is_numeric_dtype(series):
        return series.astype(float).to_numpy()
    return pd.Categorical(series).codes.astype(float)


def _load_rds(path: str) -> pd.DataFrame:
    data = pyreadr.read_r(path)
    if not data:
        raise ValueError(f"No objects found in {path}")
    obj = next(iter(data.values()))
    if not isinstance(obj, pd.DataFrame):
        raise TypeError(f"Expected a data.frame in {path}, got {type(obj)}")
    return obj


def _sample_test_df(df: pd.DataFrame, batch_col: str, sex_col: str, test_n: int, seed: int) -> pd.DataFrame:
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
    parser = argparse.ArgumentParser(
        description="Run ComBat-GAM using neuroHarmonize (native smooth_terms) on SC tables."
    )
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
    parser.add_argument("--eb", action="store_true", default=True)
    parser.add_argument("--no-eb", dest="eb", action="store_false")
    parser.add_argument("--smooth-df", type=int, default=3)
    parser.add_argument("--smooth-degree", type=int, default=3)
    parser.add_argument("--smooth-fx", action="store_true", default=True)
    parser.add_argument("--smooth-kfold", dest="smooth_fx", action="store_false")
    parser.add_argument("--ref-batch", default=None)
    args = parser.parse_args()

    df = _load_rds(args.input_rds)
    sc_cols = [c for c in df.columns if c.startswith("SC.")]
    if not sc_cols:
        raise ValueError("No SC.* columns found in input data.")

    extra_cols = [c for c in args.extra_covars.split(",") if c]
    covar_cols = [args.id_col, args.batch_col, args.age_col, args.sex_col, args.mean_fd_col] + extra_cols
    covar_cols = list(dict.fromkeys(covar_cols))

    needed_cols = sc_cols + covar_cols
    missing = [c for c in needed_cols if c not in df.columns]
    if missing:
        raise ValueError(f"Missing columns in input: {missing}")

    df = df[needed_cols].dropna()
    df = _sample_test_df(df, args.batch_col, args.sex_col, args.test_n, args.seed)

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
        eb=args.eb,
        smooth_terms=[args.age_col],
        smooth_df=args.smooth_df,
        smooth_degree=args.smooth_degree,
        smooth_fx=args.smooth_fx,
        ref_batch=args.ref_batch,
    )

    out = df[[args.id_col]].copy()
    for col in [args.batch_col, args.age_col, args.sex_col, args.mean_fd_col] + extra_cols:
        out[col] = df[col].values
    for idx, col in enumerate(sc_cols):
        out[f"{col}_h"] = data_adj[:, idx]

    out_dir = os.path.dirname(os.path.abspath(args.output_rds))
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    pyreadr.write_rds(args.output_rds, out)


if __name__ == "__main__":
    main()

