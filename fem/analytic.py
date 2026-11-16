"""
analytic.py - Beam-theory reference solutions used to VERIFY the FE results.

Cantilever (clamped at x = 0), rectangular section h x t.
  bending : pressure p on the top face   -> line load q = p * t (N/m)
  tension : pressure p on the end face   -> axial force F = p * h * t
"""


def section(h, t):
    return h * t, t * h ** 3 / 12.0          # A, I


def bending_reference(L, h, t, p, E, nu, kappa=5.0 / 6.0):
    A, I = section(h, t)
    q = p * t
    G = E / (2.0 * (1.0 + nu))
    d_bend = q * L ** 4 / (8.0 * E * I)       # Euler-Bernoulli, uniform load
    d_shear = q * L ** 2 / (2.0 * kappa * G * A)
    return {
        "tip_euler_bernoulli": d_bend,
        "tip_timoshenko": d_bend + d_shear,
        # top-fibre bending stress sigma_xx(x) = M(x) c / I,  M(x) = q (L - x)^2 / 2
        "sigma_top": lambda x: q * (L - x) ** 2 / 2.0 * (h / 2.0) / I,
        "sigma_root": q * L ** 2 / 2.0 * (h / 2.0) / I,
    }


def tension_reference(L, h, t, p, E):
    return {"elongation": p * L / E, "stress": p, "strain": p / E}
