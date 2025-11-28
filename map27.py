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


def evolve_core_coupled(phi_obj_start=0.1, phi_obj_end=4.5,
                        phi_well_start=0.5, phi_well_end=3.5,
                        steps=150):
    """
    Coupled evolution: Φ_obj and Φ_well both evolve with time.
    - Φ_obj grows (core coherence increasing)
    - Φ_well deepens (mass/flux accretion increasing well capacity)
    """
    t = np.linspace(0, 1, steps)

    # Simple linear ramps (you can swap in nonlinear later)
    phi_obj = phi_obj_start + (phi_obj_end - phi_obj_start) * t
    phi_well = phi_well_start + (phi_well_end - phi_well_start) * t

    classes = [classify_object(po, pw) for po, pw in zip(phi_obj, phi_well)]
    dphi = phi_obj - phi_well

    df = pd.DataFrame({
        "t": t,
        "Phi_obj": phi_obj,
        "Phi_well": phi_well,
        "DeltaPhi": dphi,
        "Class": classes
    })
    return df

df_coupled = evolve_core_coupled()

# Find transitions in class
transitions_coupled = df_coupled[df_coupled["Class"].ne(df_coupled["Class"].shift())]

df_coupled.head(), transitions_coupled
