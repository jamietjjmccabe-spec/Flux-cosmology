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
    phi_obj  = phi_obj_start  + (phi_obj_end  - phi_obj_start)  * t
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
