#!/usr/bin/env python3
"""
flux_priority_calibrator_census.py

Create a priority host-property census list from the Pantheon+/SH0ES mass-baseline file.

Purpose
-------
Before the full host catalogue is available, focus manual/targeted lookup on the most
important SNe:

1. SNe with duplicate/light-curve weight in Pantheon+
2. SH0ES calibrators
3. SNe used in the SH0ES Hubble-flow subset
4. Very nearby / famous low-z anchors

Outputs
-------
flux_priority_calibrator_census.csv
flux_priority_host_proxy_entry_template.csv
flux_priority_census_summary.csv

The proxy template is ready to become:

    host_structure_maturity_proxies.csv

after filling Reff_kpc, stellar_age_Gyr, sigma_v_km_s, etc.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import numpy as np
import pandas as pd


def make_priority_census(baseline_path: str | Path, outdir: str | Path, top_n: int = 150):
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(baseline_path)

    # Ensure expected flags exist.
    for col in ["IS_CALIBRATOR", "USED_IN_SH0ES_HF"]:
        if col not in df.columns:
            df[col] = 0

    # Group by unique SNID/CID. The duplicate count is useful: one host lookup applies to all duplicate entries.
    g = df.groupby("SNID", dropna=False).agg(
        CID=("CID", "first") if "CID" in df.columns else ("SNID", "first"),
        n_entries=("SNID", "size"),
        mean_z=("z_use", "mean"),
        min_z=("z_use", "min"),
        max_z=("z_use", "max"),
        logMstar=("logMstar", "mean"),
        mean_mu_resid=("mu_resid_centered", "mean"),
        std_mu_resid=("mu_resid_centered", "std"),
        mean_H0_proxy=("H0_proxy_centered", "mean"),
        is_calibrator=("IS_CALIBRATOR", "max"),
        used_in_shoes_hf=("USED_IN_SH0ES_HF", "max"),
    ).reset_index()

    # Priority score:
    # - duplicate entries matter
    # - calibrator and Hubble-flow flags matter
    # - nearby anchors matter
    # - high absolute residuals can be interesting
    g["abs_mean_mu_resid"] = g["mean_mu_resid"].abs()
    g["nearby_bonus"] = np.where(g["mean_z"] < 0.02, 1.0, 0.0)
    g["priority_score"] = (
        g["n_entries"]
        + 10.0 * g["is_calibrator"].fillna(0)
        + 5.0 * g["used_in_shoes_hf"].fillna(0)
        + 2.0 * g["nearby_bonus"]
        + 4.0 * g["abs_mean_mu_resid"].fillna(0)
    )

    g = g.sort_values(["priority_score", "n_entries", "is_calibrator", "used_in_shoes_hf"], ascending=False)

    priority = g.head(top_n).copy()

    # Host proxy entry template.
    template = priority[[
        "SNID", "CID", "n_entries", "mean_z", "logMstar",
        "is_calibrator", "used_in_shoes_hf", "mean_mu_resid", "mean_H0_proxy",
        "priority_score"
    ]].copy()

    # Add blank host fields to fill.
    template["Reff_kpc"] = ""
    template["stellar_age_Gyr"] = ""
    template["sigma_v_km_s"] = ""
    template["logsSFR"] = ""
    template["metallicity"] = ""
    template["morphology"] = ""
    template["host_name"] = ""
    template["source_radius"] = ""
    template["source_age"] = ""
    template["source_sigma"] = ""
    template["notes"] = ""

    summary = pd.DataFrame([{
        "input_rows": len(df),
        "unique_SNIDs": df["SNID"].nunique(),
        "top_n": top_n,
        "priority_rows": len(priority),
        "total_entries_covered_by_priority": int(priority["n_entries"].sum()),
        "fraction_entries_covered": float(priority["n_entries"].sum() / len(df)),
        "calibrators_in_priority": int(priority["is_calibrator"].sum()),
        "shoes_hf_in_priority": int(priority["used_in_shoes_hf"].sum()),
        "interpretation": (
            "This priority list is a partial-audit target set. Fill the blank host fields for these SNIDs first, then run flux_neighborhood_audit.py."
        ),
    }])

    priority.to_csv(outdir / "flux_priority_calibrator_census.csv", index=False)
    template.to_csv(outdir / "flux_priority_host_proxy_entry_template.csv", index=False)
    summary.to_csv(outdir / "flux_priority_census_summary.csv", index=False)

    return summary, priority, template


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--baseline", required=True)
    ap.add_argument("--outdir", default=".")
    ap.add_argument("--top-n", type=int, default=150)
    args = ap.parse_args()

    summary, priority, template = make_priority_census(args.baseline, args.outdir, args.top_n)

    print("\nFlux priority calibrator census")
    print("=" * 80)
    print(summary.to_string(index=False))
    print("\nTop 30 targets")
    cols = ["SNID", "n_entries", "mean_z", "logMstar", "is_calibrator", "used_in_shoes_hf", "mean_mu_resid", "priority_score"]
    print(priority[cols].head(30).to_string(index=False))


if __name__ == "__main__":
    main()
