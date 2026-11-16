"""
solver.py - Global assembly, boundary conditions, loads and the linear solve.

The global system is  K u = F.
Each node has 2 DOFs: u (x) at index 2*n and v (y) at index 2*n + 1.

A dense NumPy matrix is used on purpose: it keeps the code readable and runs
everywhere, including the browser. Meshes up to roughly 1,500 DOF solve in
about a second. A research code would use scipy.sparse instead.
"""

import numpy as np

from .elements import element_stiffness


def element_dofs(conn):
    """Global DOF indices of an element: [2n0, 2n0+1, 2n1, 2n1+1, ...]."""
    conn = np.asarray(conn)
    return np.column_stack([2 * conn, 2 * conn + 1]).ravel()


def assemble_stiffness(nodes, elements, D, t, n_gauss=2):
    ndof = 2 * len(nodes)
    K = np.zeros((ndof, ndof))
    for conn in elements:
        ke = element_stiffness(nodes[conn], D, t, n_gauss)
        dofs = element_dofs(conn)
        K[np.ix_(dofs, dofs)] += ke            # scatter-add into the global matrix
    return K


def edge_load_vector(nodes, edge, total_force, direction):
    """
    Consistent nodal forces for a UNIFORM traction on a straight edge whose
    resultant is total_force, acting in 'x' or 'y'.

    For linear edges each segment of length l_s carries q * l_s, shared
    half-and-half between its two end nodes.
    """
    ndof = 2 * len(nodes)
    F = np.zeros(ndof)
    pts = nodes[edge]
    seg = np.linalg.norm(np.diff(pts, axis=0), axis=1)
    q = total_force / seg.sum()               # force per unit edge length
    comp = 0 if direction == "x" else 1
    for s, length in enumerate(seg):
        for n in (edge[s], edge[s + 1]):
            F[2 * n + comp] += 0.5 * q * length
    return F


def solve(K, F, fixed_dofs, prescribed=None):
    """
    Solve K u = F with u[fixed_dofs] = prescribed (default 0).
    Returns the displacement vector u and the reaction vector R = K u - F.
    """
    ndof = K.shape[0]
    fixed_dofs = np.asarray(fixed_dofs, dtype=int)
    u = np.zeros(ndof)
    if prescribed is not None:
        u[fixed_dofs] = prescribed
    free = np.setdiff1d(np.arange(ndof), fixed_dofs)

    # partitioned system: K_ff u_f = F_f - K_fc u_c
    rhs = F[free] - K[np.ix_(free, fixed_dofs)] @ u[fixed_dofs]
    u[free] = np.linalg.solve(K[np.ix_(free, free)], rhs)

    R = K @ u - F
    return u, R
