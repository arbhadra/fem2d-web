"""
Verification tests. Run with:   python -m pytest -v

Each test checks a property that a correct finite element code MUST have.
If you change the code and a test fails, the test tells you what broke.
"""

import numpy as np
import pytest

from fem import D_matrix, load_sweep, run_case
from fem.elements import cst_stiffness, q4_stiffness
from fem.solver import assemble_stiffness, solve



# ---------------------------------------------------------------- materials
def test_D_matches_textbook():
    E, nu = 200e9, 0.3
    expected = E / (1 - nu ** 2) * np.array([[1, nu, 0], [nu, 1, 0], [0, 0, (1 - nu) / 2]])
    assert np.allclose(D_matrix(E, nu), expected)


def test_bad_material_is_rejected():
    with pytest.raises(ValueError):
        D_matrix(200e9, 0.5)


# ---------------------------------------------------------------- elements
def rigid_body_count(k, tol=1e-8):
    eig = np.linalg.eigvalsh(k)
    return int(np.sum(np.abs(eig) < tol * np.abs(eig).max()))


def D_steel():
    return D_matrix(200e9, 0.3)


def test_cst_has_three_rigid_body_modes():
    xy = np.array([[0.0, 0.0], [1.0, 0.2], [0.3, 0.9]])
    k = cst_stiffness(xy, D_steel(), 0.01)
    assert np.allclose(k, k.T)
    assert rigid_body_count(k) == 3            # 2 translations + 1 rotation


def test_q4_full_integration_has_three_rigid_body_modes():
    xy = np.array([[0.0, 0.0], [1.0, 0.1], [1.2, 1.0], [-0.1, 0.8]])
    k = q4_stiffness(xy, D_steel(), 0.01, n_gauss=2)
    assert np.allclose(k, k.T)
    assert rigid_body_count(k) == 3


def test_q4_one_point_integration_shows_hourglass_modes():
    """Reduced (1-point) integration adds 2 spurious zero-energy 'hourglass' modes."""
    xy = np.array([[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]])
    k = q4_stiffness(xy, D_steel(), 0.01, n_gauss=1)
    assert rigid_body_count(k) == 5


# ---------------------------------------------------------------- patch test
def patch_mesh(kind):
    """2x2 patch with the centre node deliberately moved off-centre."""
    nodes = np.array([[0, 0], [1, 0], [2, 0],
                      [0, 1], [1.15, 0.85], [2, 1],
                      [0, 2], [1, 2], [2, 2]], dtype=float)
    quads = np.array([[0, 1, 4, 3], [1, 2, 5, 4], [3, 4, 7, 6], [4, 5, 8, 7]])
    if kind == "Q4":
        return nodes, quads
    tris = np.vstack([quads[:, [0, 1, 2]], quads[:, [0, 2, 3]]])
    return nodes, tris


@pytest.mark.parametrize("kind", ["CST", "Q4"])
def test_patch_test_reproduces_linear_field(kind):
    """
    Impose a linear displacement field on the boundary. A convergent element
    must reproduce it EXACTLY at the interior node (constant strain state).
    """
    nodes, elems = patch_mesh(kind)
    D = D_steel()
    K = assemble_stiffness(nodes, elems, D, t=0.01)

    def exact(x, y):
        return 1e-3 * (1 + 2 * x + 3 * y), 1e-3 * (-1 + 0.5 * x - 4 * y)

    boundary = [0, 1, 2, 3, 5, 6, 7, 8]
    fixed = np.array([[2 * n, 2 * n + 1] for n in boundary]).ravel()
    values = np.array([exact(*nodes[n]) for n in boundary]).ravel()
    u, _ = solve(K, np.zeros(2 * len(nodes)), fixed, values)

    ux, uy = exact(*nodes[4])
    assert np.isclose(u[8], ux, rtol=1e-10)
    assert np.isclose(u[9], uy, rtol=1e-10)


# ---------------------------------------------------------------- full model
def test_equilibrium_reactions_balance_pressure():
    r = run_case(p=50e3, nx=20, ny=2)          # total force = p * L * t = 500 N
    assert np.isclose(r["reaction_y"], 500.0, rtol=1e-9)
    assert abs(r["reaction_x"]) < 1e-6


def test_tension_matches_pL_over_E():
    r = run_case(load="tension", p=10e6, nx=40, ny=4)
    assert abs(r["tip_error_pct"]) < 0.5      # clamped root adds a tiny Poisson effect
    # away from the clamped root the strain is uniform: eps = p / E
    assert np.isclose(np.median(r["strain"][:, 0]), 10e6 / 200e9, rtol=0.01)


def test_q4_cantilever_converges_to_beam_theory():
    r = run_case(nx=80, ny=8, element="Q4")
    assert abs(r["tip_error_pct"]) < 1.5
    assert np.isclose(r["sigma_mid"], r["sigma_mid_ref"], rtol=0.02)


def test_q4_beats_cst_on_same_mesh():
    q4 = run_case(nx=40, ny=4, element="Q4")
    cst = run_case(nx=40, ny=4, element="CST")
    assert abs(q4["tip_error_pct"]) < abs(cst["tip_error_pct"])


def test_mesh_refinement_reduces_error():
    errs = [abs(run_case(nx=n, ny=n // 10)["tip_error_pct"]) for n in (20, 40, 80)]
    assert errs[0] > errs[1] > errs[2]


def test_results_scale_linearly_with_pressure():
    r1, r2 = run_case(p=1e5), run_case(p=3e5)
    assert np.isclose(r2["tip"], 3 * r1["tip"], rtol=1e-9)
    sweep = load_sweep(r1)
    assert np.isclose(sweep["p_first_yield"] * r1["max_von_mises"] / r1["inputs"]["p"], r1["yield"])
