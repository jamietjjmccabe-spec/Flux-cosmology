import numpy as np
import pandas as pd

def classify_object(phi_obj, phi_well,
                    eps_star=0.1,
                    bd_max=0.5,
                    pulsar_min=-0.5,
                    high_phi_bh=3.0):
    """
    Tiny toy classifier for compact-object type in Φ-space.

    Parameters
    ----------
    phi_obj : float
        Object core flux-coherence (Φ_obj).
    phi_well : float
        Well flux-coherence capacity (Φ_well).
    eps_star : float
        Max |ΔΦ| for main-sequence 'matched' regime.
    bd_max : float
        Max positive ΔΦ for brown dwarf regime.
    pulsar_min : float
        Min negative ΔΦ for pulsar regime (ΔΦ < pulsar_min).
    high_phi_bh : float
        If both Φ_obj and Φ_well exceed this, classify as BH.

    Returns
    -------
    str
        One of {"BH", "pulsar/NS", "star", "brown dwarf", "subcritical"}.
    """
    dphi = phi_obj - phi_well

    # Black hole corner: both very high Φ and strongly collapsed
    if phi_obj > high_phi_bh and phi_well > high_phi_bh:
        return "BH"

    # Pulsar / NS: deep well, overpowered environment (ΔΦ strongly negative)
    if dphi < pulsar_min and phi_well > 1.0:
        return "pulsar/NS"

    # Main-sequence star: matched object and well
    if abs(dphi) <= eps_star and phi_obj >= 1.0 and phi_well >= 1.0:
        return "star"

    # Brown dwarf: demand > supply but not extreme
    if dphi > 0 and dphi <= bd_max and phi_obj >= 0.5:
        return "brown dwarf"

    # Subcritical / diffuse: below all thresholds
    return "subcritical"

# Let's make a small grid of (Φ_obj, Φ_well) and classify
phi_obj_vals = [0.3, 0.6, 1.0, 1.2, 2.0, 4.0]
phi_well_vals = [0.3, 0.6, 1.0, 1.2, 2.0, 4.0]

rows = []
for po in phi_obj_vals:
    for pw in phi_well_vals:
        label = classify_object(po, pw)
        rows.append({"Phi_obj": po, "Phi_well": pw, "DeltaPhi": po-pw, "Class": label})

df = pd.DataFrame(rows)
df_pivot = df.pivot(index="Phi_obj", columns="Phi_well", values="Class")

df, df_pivot
