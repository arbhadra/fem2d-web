"""
mesh.py - Structured meshes of a rectangle L x h.

Coordinates: x runs 0 -> L along the bar, y runs -h/2 -> +h/2 so that the
neutral axis of the beam sits at y = 0.

Node numbering is row by row:  node = j * (nx + 1) + i
    i = 0..nx  (column, along x)
    j = 0..ny  (row, along y)
"""

import numpy as np


def rectangle(L, h, nx, ny, kind="Q4"):
    """
    Return (nodes, elements).

    nodes    : (n_nodes, 2) array of x, y coordinates
    elements : (n_elem, 4) for Q4 or (n_elem, 3) for CST, node indices listed
               counter-clockwise (needed for a positive area / Jacobian)
    """
    if nx < 1 or ny < 1:
        raise ValueError("nx and ny must be at least 1")
    kind = kind.upper()

    xs = np.linspace(0.0, L, nx + 1)
    ys = np.linspace(-h / 2.0, h / 2.0, ny + 1)
    X, Y = np.meshgrid(xs, ys)                # shape (ny+1, nx+1)
    nodes = np.column_stack([X.ravel(), Y.ravel()])

    def nid(i, j):
        return j * (nx + 1) + i

    elems = []
    for j in range(ny):
        for i in range(nx):
            n0, n1 = nid(i, j), nid(i + 1, j)
            n2, n3 = nid(i + 1, j + 1), nid(i, j + 1)
            if kind == "Q4":
                elems.append([n0, n1, n2, n3])
            elif kind == "CST":
                # split each rectangle into two triangles along a diagonal
                elems.append([n0, n1, n2])
                elems.append([n0, n2, n3])
            else:
                raise ValueError("kind must be 'Q4' or 'CST'")
    return nodes, np.array(elems, dtype=int)


def edge_nodes(nodes, x=None, y=None, tol=1e-9):
    """Indices of nodes lying on the line x = const or y = const."""
    if x is not None:
        idx = np.where(np.abs(nodes[:, 0] - x) < tol)[0]
        return idx[np.argsort(nodes[idx, 1])]
    if y is not None:
        idx = np.where(np.abs(nodes[:, 1] - y) < tol)[0]
        return idx[np.argsort(nodes[idx, 0])]
    raise ValueError("give x or y")
