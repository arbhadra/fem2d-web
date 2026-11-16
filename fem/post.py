"""
post.py - Strains and stresses from displacements.

Evaluated at each element's nodes, then averaged at shared nodes
("nodal averaging") to give smooth contour plots.
"""

import numpy as np

from .elements import element_nodal_strains
from .solver import element_dofs


def nodal_fields(nodes, elements, u, D):
    """Return averaged strain (n, 3), stress (n, 3) and von Mises (n,)."""
    acc = np.zeros((len(nodes), 3))
    count = np.zeros(len(nodes))
    for conn in elements:
        acc[conn] += element_nodal_strains(nodes[conn], u[element_dofs(conn)])
        count[conn] += 1
    eps = acc / count[:, None]
    sig = eps @ D.T
    return eps, sig, von_mises(sig)


def von_mises(sig):
    """Plane-stress von Mises: sqrt(sx^2 - sx*sy + sy^2 + 3 txy^2)."""
    sx, sy, txy = sig[:, 0], sig[:, 1], sig[:, 2]
    return np.sqrt(sx ** 2 - sx * sy + sy ** 2 + 3.0 * txy ** 2)
