"""
flux_mip_v3q_sparc_pipeline.py

Parse SPARC MRT files and run v3o/v3q/v3q-retention residual tests.
Expected inputs in the working directory:
    SPARC_Lelli2016c.mrt
    MassModels_Lelli2016c.mrt or MassModels_Lelli2016c.mrt.txt
"""
from __future__ import annotations
from pathlib import Path
import io
import numpy as np
import pandas as pd
from flux_mip_v3q_model import (
    G_KPC_KMS2_PER_MSUN, SIGMA_REF, eta_v3o, memory_radius,
    remnant_proxy, retention_gate, eta_effective,
    baryonic_velocity_squared, historical_velocity_squared,
)


def get_clean_io(filename: str | Path) -> io.StringIO:
    lines = Path(filename).read_text().splitlines(True)
    start = 0
    for i, line in enumerate(lines):
        if "--------------------------------------------------------------------------------" in line:
            start = i + 1
    return io.StringIO("".join(lines[start:]))


def load_sparc_global(path="SPARC_Lelli2016c.mrt") -> pd.DataFrame:
    raw = pd.read_csv(get_clean_io(path), sep=r"\s+", header=None)
    # Columns in SPARC_Lelli2016c.mrt:
    # 0 Galaxy, 1 T, 7 L[3.6] (1e9 Lsun), 11 Rdisk kpc,
    # 12 Vflat, 13 MHI (1e8 Msun), 16 Q
    df = raw[[0, 1, 7, 11, 12, 13, 16]].copy()
    df.columns = ["Galaxy", "T", "L36", "R_disk", "Vflat", "MHI", "Q"]
    for c in df.columns[1:]:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df["M_star"] = 0.5 * df["L36"] * 1e9
    df["M_gas"] = 1.33 * df["MHI"] * 1e8
    df["M_bar"] = df["M_star"] + df["M_gas"]
    df["f_gas"] = df["M_gas"] / df["M_bar"]
    df = df[(df["R_disk"] > 0) & (df["M_bar"] > 0)].copy()
    df["Sigma_bar"] = df["M_bar"] / (2.0 * np.pi * df["R_disk"] ** 2)
    df["D_spatial"] = (df["Sigma_bar"] / SIGMA_REF) ** 0.08
    return df


def load_sparc_models(path="MassModels_Lelli2016c.mrt") -> pd.DataFrame:
    p = Path(path)
    if not p.exists() and Path(str(path) + ".txt").exists():
        p = Path(str(path) + ".txt")
    raw = pd.read_csv(get_clean_io(p), sep=r"\s+", header=None)
    # 0 Galaxy, 2 R, 3 Vobs, 4 eVobs, 5 Vgas, 6 Vdisk, 7 Vbul
    df = raw[[0, 2, 3, 4, 5, 6, 7]].copy()
    df.columns = ["Galaxy", "R", "Vobs", "eVobs", "Vgas", "Vdisk", "Vbul"]
    for c in df.columns[1:]:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    return df


def run_pipeline(global_path="SPARC_Lelli2016c.mrt", models_path="MassModels_Lelli2016c.mrt", out_dir="outputs"):
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    glob = load_sparc_global(global_path)
    mod = load_sparc_models(models_path)
    df = mod.merge(glob, on="Galaxy", how="inner")
    df = df.dropna(subset=["R", "Vobs", "Vgas", "Vdisk", "Vbul", "D_spatial", "M_bar", "R_disk", "T", "f_gas"])
    df = df[(df["R"] > 0)].copy()

    df["Vmax"] = df.groupby("Galaxy")["Vobs"].transform("max")
    df["Vbar2"] = baryonic_velocity_squared(df["Vgas"].values, df["Vdisk"].values, df["Vbul"].fillna(0).values)
    D = df["D_spatial"].values

    for label, use_ret in [("v3o", None), ("v3q", False), ("v3q_retention", True)]:
        if label == "v3o":
            eta_eff = eta_v3o(D)
        else:
            eta_eff = eta_effective(D, df["f_gas"].values, df["T"].values, vmax=df["Vmax"].values, use_retention=use_ret)
        vh2 = historical_velocity_squared(df["R"].values, df["M_bar"].values, df["R_disk"].values, D, eta_eff)
        df[f"V_{label}"] = np.sqrt(np.maximum(df["Vbar2"].values + vh2, 0.0))
        df[f"res_{label}"] = df[f"V_{label}"] - df["Vobs"]
        df[f"abs_res_{label}"] = np.abs(df[f"res_{label}"])

    galaxy = df.groupby("Galaxy").agg(
        T=("T", "first"), f_gas=("f_gas", "first"), D_spatial=("D_spatial", "first"),
        R_disk=("R_disk", "first"), M_bar=("M_bar", "first"), Vmax=("Vmax", "first"),
        MAE_v3o=("abs_res_v3o", "mean"), MAE_v3q=("abs_res_v3q", "mean"),
        MAE_v3q_retention=("abs_res_v3q_retention", "mean"),
    ).reset_index()
    summary = pd.DataFrame([
        {"model":"v3o", "point_MAE":df["abs_res_v3o"].mean(), "galaxy_MAE":galaxy["MAE_v3o"].mean()},
        {"model":"v3q", "point_MAE":df["abs_res_v3q"].mean(), "galaxy_MAE":galaxy["MAE_v3q"].mean()},
        {"model":"v3q_retention", "point_MAE":df["abs_res_v3q_retention"].mean(), "galaxy_MAE":galaxy["MAE_v3q_retention"].mean()},
    ])
    df.to_csv(out/"v3q_point_profiles.csv", index=False)
    galaxy.to_csv(out/"v3q_galaxy_summary.csv", index=False)
    summary.to_csv(out/"v3q_summary_metrics.csv", index=False)
    return df, galaxy, summary


if __name__ == "__main__":
    _, _, summary = run_pipeline()
    print(summary.to_string(index=False))