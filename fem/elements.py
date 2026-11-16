"""
elements.py - Element stiffness matrices and strain-displacement (B) matrices.

Degrees of freedom per element are ordered [u1, v1, u2, v2, ...].

CST : 3-node constant-strain triangle (linear shape functions, exact closed form)
Q4  : 4-node bilinear isoparametric quadrilateral, Gauss quadrature
"""

import numpy as np

# ---------------------------------------------------------------------------
# CST
# ---------------------------------------------------------------------------


def cst_B(xy):
    """B matrix (3 x 6) and area of a CST. xy is (3, 2), counter-clockwise."""
    (x1, y1), (x2, y2), (x3, y3) = xy
    two_A = (x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1)
    if two_A <= 0:
        raise ValueError("CST has zero or negative area (check node order).")
    b = np.array([y2 - y3, y3 - y1, y1 - y2])   # = dN_i/dx * 2A
    c = np.array([x3 - x2, x1 - x3, x2 - x1])   # = dN_i/dy * 2A
    B = np.zeros((3, 6))
    B[0, 0::2] = b
    B[1, 1::2] = c
    B[2, 0::2] = c
    B[2, 1::2] = b
    return B / two_A, two_A / 2.0


def cst_stiffness(xy, D, t):
    """k = t * A * B^T D B  (the integrand is constant, so no quadrature)."""
    B, A = cst_B(xy)
    return t * A * B.T @ D @ B


# ---------------------------------------------------------------------------
# Q4
# ---------------------------------------------------------------------------

# natural coordinates of the 4 corner nodes (counter-clockwise)
XI_NODES = np.array([[-1.0, -1.0], [1.0, -1.0], [1.0, 1.0], [-1.0, 1.0]])


def gauss_points(n):
    """1D Gauss-Legendre points and weights on [-1, 1]."""
    if n == 1:
        return np.array([0.0]), np.array([2.0])
    if n == 2:
        g = 1.0 / np.sqrt(3.0)
        return np.array([-g, g]), np.array([1.0, 1.0])
    if n == 3:
        g = np.sqrt(3.0 / 5.0)
        return np.array([-g, 0.0, g]), np.array([5.0, 8.0, 5.0]) / 9.0
    raise ValueError("n must be 1, 2 or 3")


def q4_shape(xi, eta):
    """Shape functions N (4,) and their derivatives dN/d(xi, eta) (2, 4)."""
    N = 0.25 * (1 + xi * XI_NODES[:, 0]) * (1 + eta * XI_NODES[:, 1])
    dN = np.empty((2, 4))
    dN[0] = 0.25 * XI_NODES[:, 0] * (1 + eta * XI_NODES[:, 1])   # dN/dxi
    dN[1] = 0.25 * XI_NODES[:, 1] * (1 + xi * XI_NODES[:, 0])    # dN/deta
    return N, dN


def q4_B(xy, xi, eta):
    """B matrix (3 x 8) and det(J) of a Q4 at natural point (xi, eta)."""
    _, dN_nat = q4_shape(xi, eta)
    J = dN_nat @ xy                    # 2x2 Jacobian: d(x, y)/d(xi, eta)
    detJ = np.linalg.det(J)
    if detJ <= 0:
        raise ValueError("Q4 Jacobian is not positive (distorted element or wrong node order).")
    dN = np.linalg.solve(J, dN_nat)    # dN/d(x, y), shape (2, 4)
    B = np.zeros((3, 8))
    B[0, 0::2] = dN[0]
    B[1, 1::2] = dN[1]
    B[2, 0::2] = dN[1]
    B[2, 1::2] = dN[0]
    return B, detJ


def q4_stiffness(xy, D, t, n_gauss=2):
    """k = sum_ij  w_i w_j  B^T D B  t det(J)   at the Gauss points."""
    pts, wts = gauss_points(n_gauss)
    k = np.zeros((8, 8))
    for xi, wi in zip(pts, wts):
        for eta, wj in zip(pts, wts):
            B, detJ = q4_B(xy, xi, eta)
            k += wi * wj * t * detJ * (B.T @ D @ B)
    return k


# ---------------------------------------------------------------------------
# dispatch helpers used by the assembler and post-processor
# ---------------------------------------------------------------------------


def element_stiffness(xy, D, t, n_gauss=2):
    if len(xy) == 3:
        return cst_stiffness(xy, D, t)
    return q4_stiffness(xy, D, t, n_gauss)


def element_nodal_strains(xy, u_e):
    """Strain [exx, eyy, gxy] evaluated at each node of the element."""
    if len(xy) == 3:
        B, _ = cst_B(xy)
        eps = B @ u_e
        return np.tile(eps, (3, 1))
    out = np.empty((4, 3))
    for a, (xi, eta) in enumerate(XI_NODES):
        B, _ = q4_B(xy, xi, eta)
        out[a] = B @ u_e
    return out
