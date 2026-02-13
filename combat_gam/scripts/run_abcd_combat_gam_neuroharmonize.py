#!/usr/bin/env python3
import argparse
import os
import re
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


VARIANTS = {
    "cognition": {
        "extra_covars": ["nihtbx_fluidcomp_uncorrected"],
        "baseline_only": True,
        "backfill": ["nihtbx_fluidcomp_uncorrected"],
    },
    "pfactor": {
        "extra_covars": ["GENERAL"],
        "baseline_only": False,
        "backfill": [],
    },
    "baseline_age_sex_meanfd": {
        "extra_covars": [],
        "baseline_only": True,
        "backfill": [],
    },
    "cbcl": {
        "extra_covars": ["cbcl_scr_syn_totprob_r"],
        "baseline_only": False,
        "backfill": ["cbcl_scr_syn_totprob_r"],
    },
    "comp_agecorrected_baseline": {
        "extra_covars": ["nihtbx_fluidcomp_agecorrected"],
        "baseline_only": True,
        "backfill": ["nihtbx_fluidcomp_agecorrected"],
    },
    "fluidcomp_fc_baseline": {
        "extra_covars": ["nihtbx_fluidcomp_fc"],
        "baseline_only": True,
        "backfill": ["nihtbx_fluidcomp_fc"],
    },
}


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
        sys.stderr.write(f"[WARN] pyreadr failed to read {path}: {e}\n")
        sys.stderr.write("[WARN] Falling back to rpy2 + readRDS(). This is slower but more compatible.\n")
        try:
            ro.globalenv[".rds_path"] = ro.StrVector([path])[0]
            r_obj = ro.r("as.data.frame(readRDS(.rds_path))")
            with localconverter(default_converter + pandas2ri.converter):
                df = ro.conversion.rpy2py(r_obj)
            if not isinstance(df, pd.DataFrame):
                raise ValueError(f"rpy2 readRDS did not return a DataFrame for {path}: {type(df)}")
            return df
        except Exception as e2:
            raise RuntimeError(
                f"Failed to read RDS via both pyreadr and rpy2: {path}. "
                f"pyreadr_error={e} rpy2_error={e2}"
            ) from e2


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


def _backfill_covariate(df, covar, demo_path, id_col):
    if covar in df.columns:
        return df
    if not os.path.exists(demo_path):
        raise FileNotFoundError(
            f"Missing column: {covar}. Tried to backfill from {demo_path} but file not found."
        )
    demo = pd.read_csv(demo_path)
    if id_col not in demo.columns or covar not in demo.columns:
        raise ValueError(
            f"Missing column: {covar}. Backfill failed because {demo_path} lacks {id_col} and/or {covar}."
        )
    demo = demo[[id_col, covar]].drop_duplicates(subset=[id_col])
    merged = df.merge(demo, how="left", on=id_col)
    missing_n = int(merged[covar].isna().sum())
    sys.stderr.write(
        f"[INFO] Backfilled '{covar}' from {demo_path}; rows={len(merged)}, missing_after_merge={missing_n}\n"
    )
    return merged


def _output_path(input_rds, output_dir, variant):
    base = os.path.basename(input_rds)
    base = re.sub(r"\.rds$", "", base, flags=re.IGNORECASE)
    base = re.sub(r"\.merge$", "", base, flags=re.IGNORECASE)
    return os.path.join(output_dir, f"{base}.combatgam_neuroharmonize_{variant}.rds")


def run_variant(
    df,
    sc_cols,
    variant,
    config,
    id_col,
    scan_col,
    batch_col,
    age_col,
    sex_col,
    mean_fd_col,
    event_col,
    baseline_label,
    demo_path,
    test_n,
    seed,
    output_dir,
    input_rds,
):
    for covar in config["backfill"]:
        df = _backfill_covariate(df, covar, demo_path, scan_col)

    extra_covars = config["extra_covars"]
    needed = [id_col, scan_col, batch_col, age_col, sex_col, mean_fd_col] + extra_covars
    if config["baseline_only"]:
        needed.append(event_col)
    needed = list(dict.fromkeys(sc_cols + needed))

    missing = [c for c in needed if c not in df.columns]
    if missing:
        raise ValueError(f"Missing columns for {variant}: {missing}")

    work = df[needed]
    if config["baseline_only"]:
        work = work[work[event_col] == baseline_label]
    work = work.dropna()
    work = _sample_test_df(work, batch_col, sex_col, test_n, seed)

    covars = pd.DataFrame()
    covars["SITE"] = work[batch_col].astype(str)
    covars[sex_col] = _ensure_numeric(work[sex_col])
    covars[mean_fd_col] = _ensure_numeric(work[mean_fd_col])
    for col in extra_covars:
        covars[col] = _ensure_numeric(work[col])

    basis = _build_mgcv_basis(work[age_col], age_col)
    covars = pd.concat([covars.reset_index(drop=True), basis.reset_index(drop=True)], axis=1)

    data_matrix = work[sc_cols].to_numpy()

    _, data_adj = harmonizationLearn(
        data_matrix,
        covars,
        smooth_terms=[],
    )

    out = work[[scan_col]].copy()
    out[id_col] = work[id_col]
    out[batch_col] = work[batch_col]
    out[age_col] = work[age_col]
    out[sex_col] = work[sex_col]
    out[mean_fd_col] = work[mean_fd_col]
    for col in extra_covars:
        out[col] = work[col]
    for idx, col in enumerate(sc_cols):
        out[f"{col}_h"] = data_adj[:, idx]

    os.makedirs(output_dir, exist_ok=True)
    out_path = _output_path(input_rds, output_dir, variant)
    pyreadr.write_rds(out_path, out)
    sys.stderr.write(f"[INFO] Wrote {variant} output: {out_path}\n")


def main():
    parser = argparse.ArgumentParser(
        description="Run ABCD ComBat-GAM using neuroHarmonize (HCP-D-aligned mgcv basis)."
    )
    parser.add_argument(
        "--input-rds",
        default="/ibmgpfs/cuizaixu_lab/xuxiaoyu/SC_development/interdataFolder_ABCD/"
        "SCdata_SA12_CV75_sumSCinvnode.sum.msmtcsd.merge.rds",
    )
    parser.add_argument(
        "--output-dir",
        default="outputs/results/combat_gam/abcd/baseline",
    )
    parser.add_argument("--variant", required=True, choices=sorted(VARIANTS.keys()))
    parser.add_argument("--test-n", type=int, default=0)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--demo-path", default=os.path.join("demopath", "DemodfScreenFinal.csv"))
    parser.add_argument("--id-col", default="subID")
    parser.add_argument("--scan-col", default="scanID")
    parser.add_argument("--batch-col", default="siteID")
    parser.add_argument("--age-col", default="age")
    parser.add_argument("--sex-col", default="sex")
    parser.add_argument("--mean-fd-col", default="mean_fd")
    parser.add_argument("--event-col", default="eventname")
    parser.add_argument("--baseline-label", default="baseline_year_1_arm_1")
    args = parser.parse_args()

    scdata = _load_rds(args.input_rds)
    sc_cols = [c for c in scdata.columns if c.startswith("SC.")]
    if not sc_cols:
        raise ValueError("No SC.* columns found in input data.")

    config = VARIANTS[args.variant]
    run_variant(
        scdata,
        sc_cols,
        args.variant,
        config,
        args.id_col,
        args.scan_col,
        args.batch_col,
        args.age_col,
        args.sex_col,
        args.mean_fd_col,
        args.event_col,
        args.baseline_label,
        args.demo_path,
        args.test_n,
        args.seed,
        args.output_dir,
        args.input_rds,
    )


if __name__ == "__main__":
    main()
