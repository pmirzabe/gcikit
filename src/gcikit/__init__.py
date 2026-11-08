"""gcikit: discretization-uncertainty estimates for CFD and FEA mesh studies."""
from .gci import (
    ConvergenceError,
    GCIResult,
    apparent_order,
    gci_from_cells,
    gci_three_grid,
    gci_two_grid,
    profile_gci,
    refinement_ratio,
)

__all__ = [
    "ConvergenceError",
    "GCIResult",
    "apparent_order",
    "gci_from_cells",
    "gci_three_grid",
    "gci_two_grid",
    "profile_gci",
    "refinement_ratio",
]
__version__ = "0.1.0"
