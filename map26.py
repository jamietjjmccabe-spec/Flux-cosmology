import numpy as np
import pandas as pd

def classify_object(phi_obj, phi_well,
                    eps_star=0.1,
                    bd_max=0.5,
                    pulsar_min=-0.5,
                    high_phi_bh=3.0):
    """
    Tiny toy classifier for compact-object type in Φ-space.
    Returns one of: "BH", "pulsar/NS", "star", "brown dwarf", "subcritical".
    """
    dphi = phi_obj - phi_well

    # Black hole: both very high Φ (deep, saturated well)
    if phi_obj > high_phi_bh and phi_well > high_phi_bh:
        return "BH"

    # Pulsar / NS: deep well, strongly overpowered environment (ΔΦ strongly negative)
    if dphi < pulsar_min and phi_well > 1.0:
        return "pulsar/NS"

    # Main-sequence star: matched object and well
    if abs(dphi) <= eps_star and phi_obj >= 1.0 and phi_well >= 1.0:
        return "star"

    # Brown dwarf: demand > supply but not extreme
    if dphi > 0 and dphi <= bd_max and phi_obj >= 0.5:
        return "brown dwarf"

    # Subcritical / diffuse state
    return "subcritical"


def evolve_core(phi_obj_start=0.1, phi_obj_end=4.5, phi_well=1.0, steps=100):
    """
    Simple evolution: Φ_obj grows from start to end in a fixed Φ_well.
    At each step, classify the regime with classify_object().
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
