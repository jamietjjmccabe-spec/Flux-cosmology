#!/usr/bin/env python3
"""
make_paper_figures.py

Generate publication-ready V3/V3.1 figures for:
Host Mass as a Shadow of Environmental Maturity

Expected input directory:
    paper_tables_v3_v31/

Expected files:
    table_3_v3_model_battle.csv
    table_4_v3_fixed_mass_fingerprint.csv
    table_5_v31_divergence_summary.csv
    active_100SNID_with_flux_columns.csv

Outputs:
    paper_figures_v3_v31/figure_1_model_battle.png/.pdf
    paper_figures_v3_v31/figure_2_fixed_mass_fingerprint.png/.pdf
    paper_figures_v3_v31/figure_3_distance_ladder_divergence.png/.pdf
    paper_figures_v3_v31/figure_captions.md

Run:
    python make_paper_figures.py

Optional:
    python make_paper_figures.py --tables paper_tables_v3_v31 --out paper_figures_v3_v31 --dpi 300
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


FIGSIZE_WIDE = (9.0, 5.4)
FIGSIZE_TALL = (9.0, 6.0)


def require_file(path: Path) -> Path:
    if not path.exists():
        raise FileNotFoundError(f"Required input file not found: {path}")
    return path


def save_figure(fig: plt.Figure, outdir: Path, stem: str, dpi: int = 300) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    fig.savefig(outdir / f"{stem}.png", dpi=dpi, bbox_inches="tight")
    fig.savefig(outdir / f"{stem}.pdf", bbox_inches="tight")
    plt.close(fig)


def clean_model_name(name: str) -> str:
    mapping = {
        "intercept": "Intercept\n(null)",
        "intercept_only": "Intercept\n(null)",
        "G_sat": "$G_{\\rm sat}$\n(saturation)",
        "mass_only": "Mass only\n$\\log M_*$",
        "mass_step_10": "Mass step\n$10^{10}M_\\odot$",
        "logD_reg": "$D_{\\rm reg}$",
        "capacity_state": "Capacity\nstate",
        "capacity_controls": "Capacity\ncontrols",
        "kitchen_sink": "Full structural\nmodel",
    }
    return mapping.get(str(name), str(name).replace("_", "\\n"))


def infer_model_columns(models: pd.DataFrame) -> tuple[str, str]:
    model_col_candidates = ["model", "Model", "name"]
    bic_col_candidates = ["BIC", "bic", "BIC_score"]

    model_col = next((c for c in model_col_candidates if c in models.columns), None)
    bic_col = next((c for c in bic_col_candidates if c in models.columns), None)

    if model_col is None:
        raise KeyError(f"Could not find model-name column in {models.columns.tolist()}")
    if bic_col is None:
        raise KeyError(f"Could not find BIC column in {models.columns.tolist()}")
    return model_col, bic_col


def make_figure_1_model_battle(tables_dir: Path, outdir: Path, dpi: int) -> None:
    """Figure 1: BIC model battle."""
    models_path = require_file(tables_dir / "table_3_v3_model_battle.csv")
    models = pd.read_csv(models_path)
    model_col, bic_col = infer_model_columns(models)

    # Lower BIC is better. Sort best on left.
    plot_df = models[[model_col, bic_col]].dropna().copy()
    plot_df = plot_df.sort_values(bic_col, ascending=True)
    plot_df["label"] = plot_df[model_col].map(clean_model_name)
    best_bic = plot_df[bic_col].min()
    plot_df["delta_BIC_vs_best"] = plot_df[bic_col] - best_bic

    fig, ax = plt.subplots(figsize=FIGSIZE_WIDE)
    x = np.arange(len(plot_df))
    ax.bar(x, plot_df["delta_BIC_vs_best"].to_numpy())
    ax.set_xticks(x)
    ax.set_xticklabels(plot_df["label"].tolist(), rotation=0, ha="center")
    ax.set_ylabel(r"$\Delta$BIC relative to best model")
    ax.set_title("Figure 1. V3 model battle: environmental maturity vs. mass")
    ax.axhline(0, linewidth=1)
    ax.grid(axis="y", alpha=0.25)

    # Annotate absolute BIC above bars.
    for i, (_, row) in enumerate(plot_df.iterrows()):
        ax.text(
            i,
            row["delta_BIC_vs_best"] + 0.15,
            f"BIC={row[bic_col]:.2f}",
            ha="center",
            va="bottom",
            fontsize=8,
            rotation=90 if len(plot_df) > 6 else 0,
        )

    ax.text(
        0.01,
        0.98,
        "Lower BIC is better. Bars show penalty relative to the best model.",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=9,
    )

    fig.tight_layout()
    save_figure(fig, outdir, "figure_1_model_battle", dpi=dpi)


def make_figure_2_fixed_mass_fingerprint(tables_dir: Path, outdir: Path, dpi: int) -> None:
    """Figure 2: fixed-mass fingerprint, High-F vs Low-F by mass bin."""
    fixed_path = require_file(tables_dir / "table_4_v3_fixed_mass_fingerprint.csv")
    fixed = pd.read_csv(fixed_path)

    # Normalize likely column names.
    mass_col = "mass_bin" if "mass_bin" in fixed.columns else fixed.columns[0]
    high_col = "High-F" if "High-F" in fixed.columns else None
    low_col = "Low-F" if "Low-F" in fixed.columns else None
    delta_col = "high_minus_low_H0" if "high_minus_low_H0" in fixed.columns else None

    if high_col is None or low_col is None:
        raise KeyError("table_4 must contain 'High-F' and 'Low-F' columns.")

    plot_df = fixed[[mass_col, low_col, high_col] + ([delta_col] if delta_col else [])].copy()
    plot_df = plot_df.dropna(subset=[low_col, high_col], how="all")

    x = np.arange(len(plot_df))
    width = 0.36

    fig, ax = plt.subplots(figsize=FIGSIZE_TALL)
    ax.bar(x - width / 2, plot_df[low_col], width, label="Low-F: Child/Teen")
    ax.bar(x + width / 2, plot_df[high_col], width, label="High-F: Adult/Overflowing")

    ax.set_xticks(x)
    ax.set_xticklabels(plot_df[mass_col].astype(str).tolist(), rotation=25, ha="right")
    ax.set_ylabel(r"Mean inferred $H_0$ proxy (km s$^{-1}$ Mpc$^{-1}$)")
    ax.set_xlabel(r"Host stellar-mass bin, $\log M_*$")
    ax.set_title("Figure 2. Fixed-mass maturity fingerprint")
    ax.grid(axis="y", alpha=0.25)
    ax.legend(frameon=False)

    # Draw delta annotations between pairs.
    ymin, ymax = ax.get_ylim()
    offset = (ymax - ymin) * 0.03
    for i, (_, row) in enumerate(plot_df.iterrows()):
        if pd.notna(row[low_col]) and pd.notna(row[high_col]):
            delta = row[high_col] - row[low_col]
            y = max(row[low_col], row[high_col]) + offset
            ax.plot([i - width / 2, i + width / 2], [y, y], linewidth=1)
            ax.text(i, y + offset * 0.3, f"Δ={delta:+.2f}", ha="center", va="bottom", fontsize=9)

    fig.tight_layout()
    save_figure(fig, outdir, "figure_2_fixed_mass_fingerprint", dpi=dpi)


def make_figure_3_distance_ladder_divergence(tables_dir: Path, outdir: Path, dpi: int) -> None:
    """Figure 3: distance ladder divergence across calibrator and Hubble-flow regimes."""
    div_path = require_file(tables_dir / "table_5_v31_divergence_summary.csv")
    div = pd.read_csv(div_path)

    needed = {"subset", "Adult_minus_Child_H0"}
    if not needed.issubset(div.columns):
        raise KeyError(f"table_5 must contain {needed}; found {div.columns.tolist()}")

    preferred_order = [
        "All Active",
        "Calibrators Only",
        "Hubble-flow Only",
        "HF: z < 0.025 (Local)",
        "HF: 0.025 <= z < 0.05 (Mid)",
        "HF: z >= 0.05 (High)",
        "HF: Core (|mu_resid| <= 0.2)",
        "HF: Outliers (|mu_resid| > 0.2)",
    ]

    rows = []
    for label in preferred_order:
        match = div[div["subset"].astype(str) == label]
        if not match.empty:
            rows.append(match.iloc[0])
    plot_df = pd.DataFrame(rows)

    if plot_df.empty:
        # Fallback: first 8 rows.
        plot_df = div.head(8).copy()

    labels = (
        plot_df["subset"]
        .astype(str)
        .str.replace("Hubble-flow", "HF", regex=False)
        .str.replace("Adult_minus_Child_H0", "ΔH0", regex=False)
        .tolist()
    )
    values = plot_df["Adult_minus_Child_H0"].astype(float).to_numpy()

    fig, ax = plt.subplots(figsize=FIGSIZE_TALL)
    y = np.arange(len(plot_df))
    ax.barh(y, values)
    ax.set_yticks(y)
    ax.set_yticklabels(labels)
    ax.invert_yaxis()
    ax.axvline(0, linewidth=1)
    ax.set_xlabel(r"Adult minus Child $H_0$ proxy (km s$^{-1}$ Mpc$^{-1}$)")
    ax.set_title("Figure 3. Distance-ladder divergence of the maturity signal")
    ax.grid(axis="x", alpha=0.25)

    for i, v in enumerate(values):
        x_text = v + (0.35 if v >= 0 else -0.35)
        ha = "left" if v >= 0 else "right"
        ax.text(x_text, i, f"{v:+.2f}", va="center", ha=ha, fontsize=9)

    ax.text(
        0.01,
        0.02,
        "Positive values: Adult/Overflowing hosts return higher inferred H0 proxy than Child/Filling hosts.",
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=9,
    )

    fig.tight_layout()
    save_figure(fig, outdir, "figure_3_distance_ladder_divergence", dpi=dpi)


def write_captions(outdir: Path) -> None:
    captions = """# Figure captions for Flux_Maturity_Preprint_V3

## Figure 1. V3 model battle: environmental maturity vs. mass
Bayesian Information Criterion comparison for the 100-SNID V3 audit. Lower BIC is preferred. The intercept-only model remains the absolute null winner, but the bounded saturation proxy $G_{\\rm sat}$ outperforms raw host stellar mass and the classical mass-step proxy among environmental predictors.

## Figure 2. Fixed-mass maturity fingerprint
Mean inferred $H_0$ proxy for Low-F and High-F host environments within fixed stellar-mass bins. A positive High-F minus Low-F offset indicates that structurally mature, high-registration hosts return higher inferred $H_0$ proxies than less mature hosts of comparable stellar mass, supporting the interpretation that host mass is a lower-resolution shadow of environmental maturity.

## Figure 3. Distance-ladder divergence of the maturity signal
Adult-minus-Child $H_0$ proxy offsets across the full active sample, calibrators, Hubble-flow objects, and Hubble-flow redshift/residual subgroups. The maturity signal is strongest in calibrators and local Hubble-flow objects, while it flattens or reverses in mid-redshift survey-dominated Hubble-flow objects, motivating the V3.1 resolution-limit interpretation.
"""
    (outdir / "figure_captions.md").write_text(captions, encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate V3/V3.1 paper figures.")
    parser.add_argument("--tables", default="paper_tables_v3_v31", help="Input table directory")
    parser.add_argument("--out", default="paper_figures_v3_v31", help="Output figure directory")
    parser.add_argument("--dpi", type=int, default=300, help="PNG DPI")
    args = parser.parse_args()

    tables_dir = Path(args.tables)
    outdir = Path(args.out)
    outdir.mkdir(parents=True, exist_ok=True)

    make_figure_1_model_battle(tables_dir, outdir, args.dpi)
    make_figure_2_fixed_mass_fingerprint(tables_dir, outdir, args.dpi)
    make_figure_3_distance_ladder_divergence(tables_dir, outdir, args.dpi)
    write_captions(outdir)

    print(f"Generated paper figures in: {outdir}")
    for p in sorted(outdir.iterdir()):
        print(f" - {p.name}")


if __name__ == "__main__":
    main()
