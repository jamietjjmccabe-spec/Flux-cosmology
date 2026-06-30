#!/usr/bin/env python3
"""Regenerate the June 2026 ECEE/planetary-battery candidate figures."""
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
plt.rcParams["svg.fonttype"] = "none"

HERE = Path(__file__).resolve().parent
OUT = HERE / "figures"
OUT.mkdir(exist_ok=True)

df = pd.read_csv(HERE / "nasa_ecee_usable_battery_candidates_compact.csv")
top = pd.read_csv(HERE / "nasa_ecee_usable_battery_top50_compact.csv")

plt.figure(figsize=(11, 8))
sc = plt.scatter(
    df["planet_radius_earth"], df["planet_mass_earth"],
    c=df["battery_coherence_usable_score"], s=70,
)
plt.colorbar(sc, label="battery_coherence_usable_score")
for _, row in df.nlargest(12, "battery_coherence_usable_score").iterrows():
    plt.annotate(row["name"], (row["planet_radius_earth"], row["planet_mass_earth"]), fontsize=7)
plt.xlabel("Planet radius / Earth")
plt.ylabel("Planet mass / Earth")
plt.title("Usable rocky-battery search space: mass-radius")
plt.tight_layout()
plt.savefig(OUT / "nasa_ecee_mass_radius_candidates.svg")
plt.close()

plt.figure(figsize=(11, 8))
sc = plt.scatter(
    df["insolation_earth"], df["core_to_stellar_ratio_relative"],
    c=df["usable_battery_score"], s=70,
)
plt.xscale("log")
plt.yscale("log")
plt.colorbar(sc, label="usable_battery_score")
for _, row in df.nlargest(12, "usable_battery_score").iterrows():
    plt.annotate(row["name"], (row["insolation_earth"], row["core_to_stellar_ratio_relative"]), fontsize=7)
plt.xlabel("Insolation / Earth")
plt.ylabel("Core-to-stellar ratio / Earth")
plt.title("Outer/low-flux usable-battery candidates")
plt.tight_layout()
plt.savefig(OUT / "nasa_ecee_core_to_stellar_candidates.svg")
plt.close()

plotdf = top.sort_values("battery_coherence_usable_score", ascending=True).tail(25)
plt.figure(figsize=(11, 8))
plt.barh(plotdf["name"], plotdf["battery_coherence_usable_score"])
plt.xlabel("Battery coherence usable score")
plt.title("NASA Archive ECEM/ECEE usable-battery candidates - top 25")
plt.tight_layout()
plt.savefig(OUT / "nasa_ecee_usable_battery_top25.svg")
plt.close()
