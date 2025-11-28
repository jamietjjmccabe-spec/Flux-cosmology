import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def classify_object(phi_obj, phi_well,
                    eps_star=0.1,
                    bd_max=0.5,
                    pulsar_min=-0.5,
                    high_phi_bh=3.0):
    dphi = phi_obj - phi_well

    if phi_obj > high_phi_bh and phi_well > high_phi_bh:
        return "BH"

    if dphi < pulsar_min and phi_well > 1.0:
        return "pulsar/NS"

    if abs(dphi) <= eps_star and phi_obj >= 1.0 and phi_well >= 1.0:
        return "star"

    if dphi > 0 and dphi <= bd_max and phi_obj >= 0.5:
        return "brown dwarf"

    return "subcritical"


def evolve_core(phi_obj_start=0.1, phi_obj_end=4.5, phi_well=1.0, steps=100):
    """
    Simple evolution: Φ_obj grows from start to end in a fixed Φ_well.
    At each step, classify the regime.
    """
    t = np.linspace(0, 1, steps)
    phi_obj = phi_obj_start + (phi_obj_end - phi_obj_start) * t
    phi_well_arr = np.full_like(phi_obj, phi_well, dtype=float)

    classes = [classify_object(po, phi_well) for po in phi_obj]
    dphi = phi_obj - phi_well_arr

    df = pd.DataFrame({
        "t": t,
        "Phi_obj": phi_obj,
        "Phi_well": phi_well_arr,
        "DeltaPhi": dphi,
        "Class": classes
    })
    return df

# Example: evolve a core in a modest well
df_example = evolve_core(phi_obj_start=0.1, phi_obj_end=4.5, phi_well=1.2, steps=120)

# Find transition points where class changes
transitions = df_example[df_example["Class"].ne(df_example["Class"].shift())]

df_example.head(), transitions
