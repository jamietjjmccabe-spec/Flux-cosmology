"""Sample a FITS magnification map at source RA/DEC coordinates.

Expected input catalog columns include: id, ra, dec, z or z_phot, and M_star_msun or mass/log mass.
This script is included for reproducibility of the local Phase-10 workflow.
"""
import argparse
import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy import units as u


def infer_mass_column(df):
    if "M_star_msun" in df.columns:
        return df["M_star_msun"]
    if "_logM_star" in df.columns:
        return 10 ** df["_logM_star"]
    if "mass" in df.columns:
        # Common EAZY convention: log stellar mass.
        return 10 ** df["mass"] if df["mass"].median() < 20 else df["mass"]
    raise KeyError("No mass column found: expected M_star_msun, _logM_star, or mass")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--catalog", required=True)
    ap.add_argument("--magnif-map", required=True)
    ap.add_argument("--out", default="phase10_ready_catalog.csv")
    ap.add_argument("--top-n", type=int, default=50)
    args = ap.parse_args()

    df = pd.read_csv(args.catalog)
    if "z" not in df.columns:
        df["z"] = df["z_phot"] if "z_phot" in df.columns else df["z500"]
    df["M_star_msun"] = infer_mass_column(df)
    df = df.sort_values("M_star_msun", ascending=False).head(args.top_n).copy()

    with fits.open(args.magnif_map) as hdul:
        data = hdul[0].data
        wcs = WCS(hdul[0].header)

    coords = SkyCoord(ra=df["ra"].values * u.deg, dec=df["dec"].values * u.deg)
    x, y = wcs.world_to_pixel(coords)
    mu = []
    valid = []
    for xi, yi in zip(x, y):
        ix, iy = int(round(xi)), int(round(yi))
        ok = 0 <= iy < data.shape[0] and 0 <= ix < data.shape[1]
        val = data[iy, ix] if ok else np.nan
        if not np.isfinite(val) or val <= 0:
            val = 1.0
            ok = False
        mu.append(float(val))
        valid.append(bool(ok))
    df["mu_lens"] = mu
    df["magnif_map_valid"] = valid
    df.to_csv(args.out, index=False)
    print(f"Saved {args.out}")
    print(df[["id", "z", "M_star_msun", "mu_lens"]].head(10).to_string(index=False))


if __name__ == "__main__":
    main()
