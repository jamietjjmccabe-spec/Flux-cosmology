import pandas as pd

K_CHI = 9.2
Q = 4
A_D = 12

rows = []
for regime, k in [
    ("CMB low", 1e-3),
    ("CMB high", 1e-2),
    ("RSD", 1e-1),
    ("transition", 9.2),
    ("galactic seed", 20.0),
    ("compact seed", 50.0),
]:
    G_env = (k / K_CHI) ** Q / (1.0 + (k / K_CHI) ** Q)
    rows.append({
        "regime": regime,
        "k_h_mpc": k,
        "G_env": G_env,
        "A_D_times_G_env": A_D * G_env,
        "A_eff_if_Gbound_1": 1 + A_D * G_env,
    })

df = pd.DataFrame(rows)
df.to_csv("flux_containment_audit.csv", index=False)
print(df.to_string(index=False))
