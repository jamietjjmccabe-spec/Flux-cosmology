#!/usr/bin/env python3
"""Compute one-sided Flux-MIP gamma bounds from a frozen exposure X.

For V/V_standard = exp(-gamma_MIP X), a one-sided unexplained fractional
visibility-loss limit epsilon gives

    gamma_max = -log(1-epsilon)/X.

This utility does not turn a symmetric one-sigma error bar into a one-sided
limit. The supplied loss limit must come from a declared likelihood or
confidence construction.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any


def load_x(path: Path) -> tuple[float, dict[str, Any]]:
    data = json.loads(path.read_text(encoding="utf-8"))
    keys = (
        "X_half_dephasing_exposure_s",
        "X_s",
        "X",
    )
    for key in keys:
        if key in data:
            return float(data[key]), data
    raise KeyError(f"No exposure key found; expected one of {keys}")


def gamma_upper_limit(x_s: float, fractional_loss: float) -> float:
    if x_s <= 0:
        raise ValueError("X must be positive.")
    if not 0 < fractional_loss < 1:
        raise ValueError("fractional_loss must lie strictly between 0 and 1.")
    return -math.log1p(-fractional_loss) / x_s


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--results", type=Path, help="JSON containing X.")
    source.add_argument("--x", type=float, help="Exposure X in seconds.")
    parser.add_argument("--max-fractional-loss", type=float, required=True)
    parser.add_argument("--output", type=Path, default=Path("flux_mip_bound.json"))
    args = parser.parse_args()

    if args.results:
        x_s, source_data = load_x(args.results)
        provenance: dict[str, Any] = {"results_file": str(args.results), "source": source_data}
    else:
        x_s = float(args.x)
        provenance = {"direct_x_input_s": x_s}

    gamma_max = gamma_upper_limit(x_s, args.max_fractional_loss)
    output = {
        "model": "V/V_standard = exp(-gamma_MIP X)",
        "X_half_dephasing_exposure_s": x_s,
        "one_sided_max_fractional_loss": args.max_fractional_loss,
        "gamma_mip_upper_limit_s_inv": gamma_max,
        "provenance": provenance,
        "warning": "The loss input must be a one-sided residual-contrast limit, not a raw symmetric error bar.",
    }
    args.output.write_text(json.dumps(output, indent=2), encoding="utf-8")
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
