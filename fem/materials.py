"""
materials.py - Material library and the plane-stress D matrix.

Isotropic, linear elastic. stress = D @ strain with
    stress = [sigma_xx, sigma_yy, tau_xy], strain = [eps_xx, eps_yy, gamma_xy]
"""

import numpy as np

# Typical handbook values (room temperature). E in Pa, yield in Pa.
# Printed polymers: tensile strength used as "yield"; they depend on print settings.
PRESETS = {
    "steel_a36":    {"label": "Structural steel (A36)",      "E": 200e9, "nu": 0.26, "yield": 250e6},
    "ss304":        {"label": "Stainless steel 304",         "E": 193e9, "nu": 0.29, "yield": 215e6},
    "al6061":       {"label": "Aluminium 6061-T6",           "E": 68.9e9, "nu": 0.33, "yield": 276e6},
    "ti64":         {"label": "Titanium Ti-6Al-4V",          "E": 113.8e9, "nu": 0.34, "yield": 880e6},
    "pla":          {"label": "PLA (3D printed, typical)",   "E": 3.5e9, "nu": 0.36, "yield": 50e6},
    "abs":          {"label": "ABS (3D printed, typical)",   "E": 2.1e9, "nu": 0.35, "yield": 40e6},
}


def D_matrix(E, nu):
    """Plane-stress constitutive matrix."""
    if E <= 0 or not (-1.0 < nu < 0.5):
        raise ValueError("Need E > 0 and -1 < nu < 0.5")
    return E / (1.0 - nu ** 2) * np.array([[1.0, nu, 0.0],
                                           [nu, 1.0, 0.0],
                                           [0.0, 0.0, (1.0 - nu) / 2.0]])
