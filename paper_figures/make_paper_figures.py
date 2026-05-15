#!/usr/bin/env python3
"""
make_paper_figures.py

Generate publication figures for the Flux/MIP V3 / V3.1 maturity-step preprint.

Inputs expected in the same directory as this script, or supplied via --indir:
  - paper_tables_v3_v31_table_3_v3_model_battle.csv
  - table_4_v3_fixed_mass_fingerprint.csv
  - paper_tables_v3_v31_table_5_v31_divergence_summary.csv
  - paper_tables_v3_v31_table_6_v31_redshift_split.csv
  - paper_tables_v3_v31_table_7_v31_calibrator_vs_hf_mass_bins.csv

Outputs:
  - figure_1_model_battle_bic.png / .pdf
  - figure_2_fixed_mass_fingerprint.png / .pdf
  - figure_3_distance_ladder_divergence.png / .pdf
"""

from pathlib import Path
import argparse
import textwrap
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def wrap_labels(labels, width=18):
    return ["\n".join(textwrap.wrap(str(x), width=width)) for x in labels]


def save_figure(fig, outpath_base):
    fig.tight_layout()
    fig.savefig(str(outpath_base) + ".png", dpi=300, bbox_inches="tight")
    fig.savefig(str(outpath_base) + ".pdf", bbox_inches="tight")
    plt.close(fig)


def figure_1_model_battle(indir: Path, outdir: Path):
    df = pd.read_csv(indir / "paper_tables_v3_v31_table_3_v3_model_battle.csv").copy()
    df = df.sort_values("BIC", ascending=True)

    fig, ax = plt.subplots(figsize=(9, 5.5))
    ax.bar(df["model"], df["BIC"])
    ax.set_title("Figure 1. V3 Model Battle: BIC Ranking")
    ax.set_ylabel("BIC (lower is better)")
    ax.set_xlabel("Model")
    ax.set_xticklabels(wrap_labels(df["model"], 12), rotation=0)

    for idx, val in enumerate(df["BIC"]):
        ax.text(idx, val, f"{val:.2f}", ha="center", va="bottom", fontsize=8)

    note = "Null remains absolute BIC winner; G_sat is the strongest environmental predictor and outperforms mass-only."
    ax.text(0.5, -0.24, note, ha="center", va="top", transform=ax.transAxes, fontsize=9)

    save_figure(fig, outdir / "figure_1_model_battle_bic")


def figure_2_fixed_mass_fingerprint(indir: Path, outdir: Path):
    df = pd.read_csv(indir / "table_4_v3_fixed_mass_fingerprint.csv").copy()

    x = np.arange(len(df))
    width = 0.36

    fig, ax = plt.subplots(figsize=(10, 5.8))
    ax.bar(x - width/2, df["Low-F"], width, label="Low-F / Child+Teen")
    ax.bar(x + width/2, df["High-F"], width, label="High-F / Adult")

    ax.set_title("Figure 2. Fixed-Mass Fingerprint")
    ax.set_ylabel("Mean inferred H0 proxy")
    ax.set_xlabel("Host stellar-mass bin")
    ax.set_xticks(x)
    ax.set_xticklabels(wrap_labels(df["mass_bin"], 14), rotation=0)
    ax.legend()

    for idx, row in df.iterrows():
        delta = row.get("high_minus_low_H0", np.nan)
        if pd.notna(delta):
            y = np.nanmax([row.get("High-F", np.nan), row.get("Low-F", np.nan)])
            ax.text(idx, y + 0.6, f"Δ={delta:+.2f}", ha="center", fontsize=8)

    note = "Within fixed mass bins, High-F hosts generally return higher H0 proxies than Low-F hosts."
    ax.text(0.5, -0.23, note, ha="center", va="top", transform=ax.transAxes, fontsize=9)

    save_figure(fig, outdir / "figure_2_fixed_mass_fingerprint")


def figure_3_distance_ladder_divergence(indir: Path, outdir: Path):
    df = pd.read_csv(indir / "paper_tables_v3_v31_table_5_v31_divergence_summary.csv").copy()

    keep = [
        "All Active",
        "Calibrators Only",
        "Hubble-flow Only",
        "HF: z < 0.025 (Local)",
        "HF: 0.025 <= z < 0.05 (Mid)",
        "HF: Core (|mu_resid| <= 0.2)",
        "HF: Outliers (|mu_resid| > 0.2)",
    ]
    plot_df = df[df["subset"].isin(keep)].copy()
    plot_df["subset"] = pd.Categorical(plot_df["subset"], categories=keep, ordered=True)
    plot_df = plot_df.sort_values("subset")

    fig, ax = plt.subplots(figsize=(11, 6.2))
    x = np.arange(len(plot_df))
    ax.bar(x, plot_df["Adult_minus_Child_H0"])
    ax.axhline(0, linewidth=1)

    ax.set_title("Figure 3. Distance-Ladder Divergence")
    ax.set_ylabel("Adult-minus-Child H0 proxy")
    ax.set_xlabel("Subset")
    ax.set_xticks(x)
    ax.set_xticklabels(wrap_labels(plot_df["subset"], 18), rotation=0)

    for idx, val in enumerate(plot_df["Adult_minus_Child_H0"]):
        if pd.notna(val):
            va = "bottom" if val >= 0 else "top"
            ax.text(idx, val, f"{val:+.2f}", ha="center", va=va, fontsize=8)

    note = "Maturity signal is strong in calibrators and local Hubble-flow, but flattens/reverses in mid-z survey flow."
    ax.text(0.5, -0.27, note, ha="center", va="top", transform=ax.transAxes, fontsize=9)

    save_figure(fig, outdir / "figure_3_distance_ladder_divergence")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--indir", default=".", help="Directory containing paper table CSV files")
    parser.add_argument("--outdir", default="publication_figures", help="Output directory for figures")
    args = parser.parse_args()

    indir = Path(args.indir)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    figure_1_model_battle(indir, outdir)
    figure_2_fixed_mass_fingerprint(indir, outdir)
    figure_3_distance_ladder_divergence(indir, outdir)

    print(f"Figures written to: {outdir}")


if __name__ == "__main__":
    main()
