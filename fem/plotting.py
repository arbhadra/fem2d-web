"""
plotting.py - Matplotlib figures used both locally and on the website.
"""

import matplotlib

matplotlib.use("Agg")                         # no window needed; works in the browser too
import matplotlib.pyplot as plt               # noqa: E402
import matplotlib.tri as mtri                 # noqa: E402
import numpy as np                            # noqa: E402
from matplotlib.collections import PolyCollection  # noqa: E402

FIELDS = {   # key: (label, getter, scale to display units)
    "von_mises": ("von Mises stress (MPa)", lambda r: r["von_mises"], 1e-6),
    "sxx": ("Bending/axial stress σxx (MPa)", lambda r: r["stress"][:, 0], 1e-6),
    "txy": ("Shear stress τxy (MPa)", lambda r: r["stress"][:, 2], 1e-6),
    "exx": ("Strain εxx (micro-strain)", lambda r: r["strain"][:, 0], 1e6),
}


def _triangles(elements):
    """Split Q4 elements into 2 triangles each so matplotlib can contour them."""
    if elements.shape[1] == 3:
        return elements
    return np.vstack([elements[:, [0, 1, 2]], elements[:, [0, 2, 3]]])


def auto_scale(res, fraction=0.08):
    """Scale factor so the largest displacement is drawn as ~8% of the bar length."""
    umax = np.abs(res["u"]).max()
    return 1.0 if umax == 0 else fraction * res["inputs"]["L"] / umax


def plot_field(res, field="von_mises"):
    """Contour of a nodal field on the DEFORMED shape; undeformed mesh in grey."""
    label, get, k = FIELDS[field]
    nodes, elems = res["nodes"], res["elements"]
    scale = auto_scale(res)
    xy = nodes + scale * res["u"].reshape(-1, 2)
    values = get(res) * k

    fig, ax = plt.subplots(figsize=(9, 3))
    ax.add_collection(PolyCollection(nodes[elems], facecolor="none",
                                     edgecolor="#bbbbbb", linewidth=0.5))
    tri = mtri.Triangulation(xy[:, 0], xy[:, 1], _triangles(elems))
    if field == "von_mises":
        tpc = ax.tripcolor(tri, values, shading="gouraud", cmap="viridis")
    else:                                     # signed field: symmetric colour scale
        vmax = np.abs(values).max() or 1.0
        tpc = ax.tripcolor(tri, values, shading="gouraud", cmap="RdBu_r", vmin=-vmax, vmax=vmax)
    fig.colorbar(tpc, ax=ax, label=label, shrink=0.85)
    ax.set_title(f"{label} on deformed shape (displacements × {scale:,.0f})")
    ax.set_aspect("equal")
    ax.autoscale()
    ax.set_xlabel("x (m)")
    ax.set_ylabel("y (m)")
    fig.tight_layout()
    return fig


def plot_top_stress(res):
    """Bending stress along the top surface: FEM vs beam theory M c / I."""
    fig, ax = plt.subplots(figsize=(6.5, 4))
    ax.plot(res["top_x"], res["top_sigma"] / 1e6, "o", ms=3, label="FEM")
    ax.plot(res["top_x"], res["top_sigma_ref"] / 1e6, "-", label="Beam theory  M·c / I")
    ax.set_xlabel("Distance from clamped end x (m)")
    ax.set_ylabel("σxx on top surface (MPa)")
    ax.set_title("Stress along the top surface")
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    return fig


def plot_load_sweep(sweep, yield_strength):
    """Applied pressure vs max von Mises stress, with the yield limit."""
    fig, ax = plt.subplots(figsize=(6.5, 4))
    p_max = max(sweep["p"][-1], sweep["p_first_yield"] * 1.1)
    ps = np.linspace(0, p_max, 50)
    slope = sweep["max_vm"][-1] / sweep["p"][-1]
    ax.plot(ps / 1e3, slope * ps / 1e6, "-", label="Max von Mises stress (FEM)")
    ax.axhline(yield_strength / 1e6, color="#d9534f", ls="--", label="Yield strength")
    ax.axvline(sweep["p_first_yield"] / 1e3, color="#888", ls=":")
    ax.plot(sweep["p"][-1] / 1e3, sweep["max_vm"][-1] / 1e6, "ko", label="Your load")
    ax.annotate(f"first yield at {sweep['p_first_yield']/1e3:,.0f} kPa",
                (sweep["p_first_yield"] / 1e3, yield_strength / 1e6),
                textcoords="offset points", xytext=(-150, 12))
    ax.set_xlabel("Applied pressure p (kPa)")
    ax.set_ylabel("Stress (MPa)")
    ax.set_title("Pressure vs stress: how much load before yielding")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="upper left")
    fig.tight_layout()
    return fig


def plot_convergence(study):
    from matplotlib.ticker import NullFormatter, ScalarFormatter
    fig, ax = plt.subplots(figsize=(6.5, 4.2))
    for el, marker in (("CST", "s-"), ("Q4", "o-")):
        d = study[el]
        ax.loglog(d[:, 0], np.abs(d[:, 3]), marker, label=el)
    for axis in (ax.xaxis, ax.yaxis):
        axis.set_major_formatter(ScalarFormatter())
        axis.set_minor_formatter(NullFormatter())
    ax.set_xlabel("Number of degrees of freedom")
    ax.set_ylabel("|error| in tip deflection vs beam theory (%)")
    ax.set_title("Mesh convergence: Q4 reaches 1% far sooner than CST")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend()
    fig.tight_layout()
    return fig
