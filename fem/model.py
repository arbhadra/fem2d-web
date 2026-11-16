"""
model.py - run_case(): one complete analysis. The website calls only this.
"""

import time

import numpy as np

from . import analytic, mesh, post, solver
from .materials import D_matrix


def run_case(L=1.0, h=0.1, t=0.01, p=100e3, load="bending",
             E=200e9, nu=0.26, yield_strength=250e6,
             nx=40, ny=4, element="Q4"):
    """
    Cantilever bar, left edge (x = 0) clamped.
    load = "bending": pressure p (Pa) pushes down on the top face
    load = "tension": pressure p (Pa) pulls on the right end face
    """
    for name, val in (("L", L), ("h", h), ("t", t)):
        if val <= 0:
            raise ValueError(f"{name} must be positive")

    t0 = time.perf_counter()
    nodes, elems = mesh.rectangle(L, h, nx, ny, element)
    D = D_matrix(E, nu)
    K = solver.assemble_stiffness(nodes, elems, D, t)

    left = mesh.edge_nodes(nodes, x=0.0)
    top = mesh.edge_nodes(nodes, y=h / 2)
    right = mesh.edge_nodes(nodes, x=L)
    fixed = np.concatenate([2 * left, 2 * left + 1])

    if load == "bending":
        F = solver.edge_load_vector(nodes, top, -p * L * t, "y")
    elif load == "tension":
        F = solver.edge_load_vector(nodes, right, p * h * t, "x")
    else:
        raise ValueError("load must be 'bending' or 'tension'")

    u, R = solver.solve(K, F, fixed)
    eps, sig, vm = post.nodal_fields(nodes, elems, u, D)
    U = u.reshape(-1, 2)

    res = {
        "nodes": nodes, "elements": elems, "u": u,
        "strain": eps, "stress": sig, "von_mises": vm,
        "n_nodes": len(nodes), "n_elements": len(elems), "n_dof": 2 * len(nodes),
        "reaction_x": R[0::2][left].sum(), "reaction_y": R[1::2][left].sum(),
        "max_von_mises": vm.max(), "max_strain": np.abs(eps[:, 0]).max(),
        "safety_factor": yield_strength / vm.max(),
        "yield": yield_strength, "solve_time_s": time.perf_counter() - t0,
        "inputs": dict(L=L, h=h, t=t, p=p, load=load, E=E, nu=nu,
                       nx=nx, ny=ny, element=element),
    }

    if load == "bending":
        ref = analytic.bending_reference(L, h, t, p, E, nu)
        tip = -U[right, 1].mean()
        res.update(tip=tip, tip_ref=ref["tip_timoshenko"],
                   top_x=nodes[top, 0], top_sigma=sig[top, 0],
                   top_sigma_ref=ref["sigma_top"](nodes[top, 0]))
        mid = top[np.argmin(np.abs(nodes[top, 0] - L / 2))]
        res["sigma_mid"] = sig[mid, 0]
        res["sigma_mid_ref"] = ref["sigma_top"](L / 2)
    else:
        ref = analytic.tension_reference(L, h, t, p, E)
        tip = U[right, 0].mean()
        res.update(tip=tip, tip_ref=ref["elongation"])
    res["tip_error_pct"] = 100 * (res["tip"] - res["tip_ref"]) / res["tip_ref"]
    return res


def load_sweep(res, n=11):
    """
    Linear elasticity: every result scales with p. So one solve gives the whole
    pressure-response curve, and first yield is at p * (yield / max von Mises).
    """
    p = res["inputs"]["p"]
    ps = np.linspace(0.0, p, n)
    return {"p": ps, "tip": res["tip"] * ps / p,
            "max_vm": res["max_von_mises"] * ps / p,
            "p_first_yield": p * res["safety_factor"]}


def convergence_study(meshes=((10, 1), (20, 2), (40, 4), (80, 8)), **kwargs):
    """Tip-deflection error vs DOF for CST and Q4 (used for verification)."""
    out = {}
    for el in ("CST", "Q4"):
        rows = []
        for nx, ny in meshes:
            r = run_case(nx=nx, ny=ny, element=el, **kwargs)
            rows.append((r["n_dof"], r["tip"], r["tip_ref"], r["tip_error_pct"]))
        out[el] = np.array(rows)
    return out
