"""fem - a small, readable 2D plane-stress finite element solver (CST and Q4)."""

from .materials import PRESETS, D_matrix  # noqa: F401
from .model import convergence_study, load_sweep, run_case  # noqa: F401
