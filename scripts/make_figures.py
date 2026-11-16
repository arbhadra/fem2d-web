"""
make_figures.py - Run the default case locally and save PNGs into assets/.

    python scripts/make_figures.py
"""

import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from fem import convergence_study, load_sweep, run_case  # noqa: E402
from fem.plotting import (plot_convergence, plot_field, plot_load_sweep,  # noqa: E402
                          plot_top_stress)

OUT = os.path.join(os.path.dirname(__file__), "..", "assets")
os.makedirs(OUT, exist_ok=True)


def save(fig, name):
    fig.savefig(os.path.join(OUT, name), dpi=150)
    print("saved assets/" + name)


if __name__ == "__main__":
    r = run_case(nx=80, ny=8)
    print(f"Steel, Q4 80x8: tip = {r['tip']*1e3:.4f} mm, beam theory = {r['tip_ref']*1e3:.4f} mm, "
          f"error = {r['tip_error_pct']:.2f} %, safety factor = {r['safety_factor']:.2f}")
    for f in ("von_mises", "sxx", "exx"):
        save(plot_field(r, f), f"{f}.png")
    save(plot_top_stress(r), "top_stress.png")
    save(plot_load_sweep(load_sweep(r), r["yield"]), "load_sweep.png")
    save(plot_convergence(convergence_study()), "convergence.png")
