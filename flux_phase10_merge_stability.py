import math
import numpy as np
import pandas as pd

files = {
    "06": "phase10_ready_catalog_z06.csv",
    "08": "phase10_ready_catalog_z08.csv",
    "10": "phase10_ready_catalog_z10.csv",
    "15": "phase10_ready_catalog_z15.csv",
    "20": "phase10_ready_catalog_z20.csv",
}

EPSILON = 0.20
A_D = 12.0
F_GATE = 1e-4
P_GATE = 1.0


def erfc_vec(x):
    try:
        from scipy.special import erfc
        return erfc(x)
    except Exception:
        return np.vectorize(math.erfc)(x)


def sigma0_proxy(M):
    logM = np.log10(np.maximum(np.asarray(M, dtype=float), 1e-30))
    x = np.array([5, 6, 7, 8, 9, 10], dtype=float)
    y = np.array([8.399696, 7.340263, 6.319449, 5.343660, 4.419663, 3.555600], dtype=float)
    out = np.interp(logM, x, y)
    high = logM > 10
    low = logM < 5
    out[high] = y[-1] + (y[-1] - y[-2]) * (logM[high] - 10)
    out[low] = y[0] + (y[1] - y[0]) * (logM[low] - 5)
    return np.maximum(out, 0.05)


def D_proxy(z):
    return 0.98 / np.power(1.0 + np.asarray(z, dtype=float), 0.86)


def G_bound_from_gate_mass(M_gate, z, f_gate=F_GATE, p=P_GATE):
    delta_c = 1.686
    sigma_z = sigma0_proxy(M_gate) * D_proxy(z)
    f_coll = erfc_vec(delta_c / (np.sqrt(2.0) * np.maximum(sigma_z, 1e-12)))
    G = np.power(f_coll / (f_coll + f_gate), p)
    return np.clip(G, 0, 1)


def classify_near_planes(row):
    z_obs = float(row["z"])
    valid = []
    for key in ["06", "08", "10", "15", "20"]:
        z_plane = int(key)
        g = row.get(f"G_{key}", np.nan)
        mu = row.get(f"mu_{key}", np.nan)
        if pd.notna(g) and pd.notna(mu) and mu > 1.0 and abs(z_plane - z_obs) <= 5.0:
            valid.append(float(g))
    if len(valid) == 0:
        return "outside_map_or_invalid"
    if all(g > 0.5 for g in valid):
        return "stable_registered"
    if any(g > 0.5 for g in valid) and any(g <= 0.5 for g in valid):
        return "caustic_sensitive_registered"
    if all(0.1 < g <= 0.5 for g in valid):
        return "weak_registered"
    if max(valid) > 0.1:
        return "weak_or_mixed"
    return "unregistered"


def main():
    dfs = {}
    for key, path in files.items():
        df = pd.read_csv(path)
        required = ["id", "z", "M_star_msun", "mu_lens"]
        missing = [c for c in required if c not in df.columns]
        if missing:
            raise KeyError(f"{path} missing columns: {missing}")
        df = df[["id", "z", "M_star_msun", "mu_lens"]].copy()
        dfs[key] = df.rename(columns={"mu_lens": f"mu_{key}"})

    master = dfs["06"]
    for key in ["08", "10", "15", "20"]:
        master = master.merge(dfs[key][["id", f"mu_{key}"]], on="id", how="left")

    for key in ["06", "08", "10", "15", "20"]:
        master[f"M_true_{key}"] = master["M_star_msun"] / master[f"mu_{key}"]
        master[f"M_gate_{key}"] = master[f"M_true_{key}"] / EPSILON
        master[f"G_{key}"] = G_bound_from_gate_mass(master[f"M_gate_{key}"], master["z"])
        master[f"flux_adv_{key}"] = 1.0 + A_D * master[f"G_{key}"]

    g_cols = [f"G_{k}" for k in ["06", "08", "10", "15", "20"]]
    master["G_min_all"] = master[g_cols].min(axis=1)
    master["G_max_all"] = master[g_cols].max(axis=1)
    master["stability_class"] = master.apply(classify_near_planes, axis=1)
    master.to_csv("flux_uncover_final_stability.csv", index=False)
    print(master["stability_class"].value_counts())


if __name__ == "__main__":
    main()
