"""Grid Convergence Index (GCI) calculations for mesh-refinement studies.

Implements the three-grid procedure in

    I. B. Celik, U. Ghia, P. J. Roache, C. J. Freitas, H. Coleman, P. E. Raad,
    "Procedure for Estimation and Reporting of Uncertainty Due to
    Discretization in CFD Applications", J. Fluids Eng. 130(7), 078001, 2008.

and the two-grid GCI of Roache (1994, 1998) with an assumed order.

Grid 1 is always the finest grid and grid 3 the coarsest.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import List, Optional, Sequence


class ConvergenceError(ValueError):
    """Raised when a study cannot give a meaningful GCI."""


def refinement_ratio(n_fine: float, n_coarse: float, dim: int) -> float:
    """Representative refinement ratio h_coarse / h_fine from cell counts.

    Uses h = (1/N)^(1/dim), so r = (N_fine / N_coarse)^(1/dim).
    """
    if dim not in (1, 2, 3):
        raise ValueError("dim must be 1, 2 or 3")
    if n_fine <= 0 or n_coarse <= 0:
        raise ValueError("cell counts must be positive")
    if n_fine <= n_coarse:
        raise ValueError("the fine grid must have more cells than the coarse grid")
    return (n_fine / n_coarse) ** (1.0 / dim)


@dataclass
class GCIResult:
    """Outcome of a three-grid study (Celik et al. 2008 notation)."""

    phi: tuple            # (phi1, phi2, phi3), fine to coarse
    r21: float
    r32: float
    p: float              # apparent order
    s: int                # sign of eps32/eps21; -1 means oscillatory convergence
    phi_ext21: float      # Richardson-extrapolated value
    phi_ext32: float
    e_a21: float          # approximate relative error, fine pair
    e_a32: float
    e_ext21: float        # extrapolated relative error, fine grid
    gci_fine21: float     # fine-grid convergence index (fraction, not %)
    gci_medium32: float
    asymptotic_ratio: Optional[float]  # GCI32 / (r21^pf * GCI21) with the formal order pf; ~1 in the asymptotic range
    iterations: int
    warnings: List[str] = field(default_factory=list)

    @property
    def oscillatory(self) -> bool:
        return self.s < 0

    def summary(self) -> str:
        lines = [
            f"phi (fine->coarse) : {self.phi[0]:.6g}, {self.phi[1]:.6g}, {self.phi[2]:.6g}",
            f"r21, r32           : {self.r21:.4f}, {self.r32:.4f}",
            f"apparent order p   : {self.p:.4f}" + ("  (oscillatory)" if self.oscillatory else ""),
            f"phi_ext21          : {self.phi_ext21:.6g}",
            f"e_a21              : {100 * self.e_a21:.3f} %",
            f"e_ext21            : {100 * self.e_ext21:.3f} %",
            f"GCI_fine21         : {100 * self.gci_fine21:.3f} %",
            f"GCI_medium32       : {100 * self.gci_medium32:.3f} %",
        ]
        if self.asymptotic_ratio is not None:
            lines.append(f"asymptotic ratio   : {self.asymptotic_ratio:.4f}")
        lines += [f"warning: {w}" for w in self.warnings]
        return "\n".join(lines)


def apparent_order(eps21: float, eps32: float, r21: float, r32: float,
                   tol: float = 1e-10, max_iter: int = 200):
    """Solve Celik eqs. (3a-c) for the apparent order p by fixed-point iteration.

    Returns (p, s, iterations).
    """
    if eps21 == 0.0 or eps32 == 0.0:
        raise ConvergenceError(
            "two grids give identical values; the apparent order is undefined")
    ratio = eps32 / eps21
    s = 1 if ratio > 0 else -1
    ln_ratio = math.log(abs(ratio))
    ln_r21 = math.log(r21)

    p = abs(ln_ratio) / ln_r21          # start from q = 0
    for it in range(1, max_iter + 1):
        num = r21 ** p - s
        den = r32 ** p - s
        if num <= 0 or den <= 0:
            raise ConvergenceError(
                "q(p) is undefined for this data (r^p - s <= 0); check the study")
        q = math.log(num / den)
        p_new = abs(ln_ratio + q) / ln_r21
        if abs(p_new - p) < tol:
            return p_new, s, it
        p = p_new
    raise ConvergenceError(f"apparent order did not converge in {max_iter} iterations")


def gci_three_grid(phi1: float, phi2: float, phi3: float,
                   r21: float, r32: float, safety_factor: float = 1.25,
                   formal_order: Optional[float] = None) -> GCIResult:
    """Three-grid GCI per Celik et al. (2008).

    phi1 is the solution on the finest grid. r21 = h2/h1 and r32 = h3/h2,
    both greater than 1. If the formal order of the scheme is given, the
    asymptotic-range ratio GCI32 / (r21^pf GCI21) is computed with it (with the
    observed order the ratio is close to 1 by construction, so it tells nothing).
    """
    if r21 <= 1.0 or r32 <= 1.0:
        raise ValueError("refinement ratios must be greater than 1 (h_coarse / h_fine)")

    eps21 = phi2 - phi1
    eps32 = phi3 - phi2
    p, s, iters = apparent_order(eps21, eps32, r21, r32)

    a21 = r21 ** p
    a32 = r32 ** p
    phi_ext21 = (a21 * phi1 - phi2) / (a21 - 1.0)
    phi_ext32 = (a32 * phi2 - phi3) / (a32 - 1.0)

    warnings: List[str] = []
    if phi1 == 0.0 or phi2 == 0.0 or phi_ext21 == 0.0:
        raise ConvergenceError("relative errors are undefined when a value is zero")

    e_a21 = abs((phi1 - phi2) / phi1)
    e_a32 = abs((phi2 - phi3) / phi2)
    e_ext21 = abs((phi_ext21 - phi1) / phi_ext21)
    gci21 = safety_factor * e_a21 / (a21 - 1.0)
    gci32 = safety_factor * e_a32 / (a32 - 1.0)
    asym: Optional[float] = None
    if formal_order is not None:
        if formal_order <= 0:
            raise ValueError("formal_order must be positive")
        f21 = r21 ** formal_order
        f32 = r32 ** formal_order
        asym = (e_a32 / (f32 - 1.0)) / (f21 * e_a21 / (f21 - 1.0))

    if s < 0:
        warnings.append("oscillatory convergence (eps32/eps21 < 0)")
    if r21 < 1.3 or r32 < 1.3:
        warnings.append("refinement ratio below 1.3; Celik et al. recommend r > 1.3")
    if not 0.5 <= p <= 6.0:
        warnings.append(f"apparent order {p:.2f} is far from typical scheme orders")
    if asym is not None and abs(asym - 1.0) > 0.1:
        warnings.append("grids may not be in the asymptotic range "
                        "(GCI32 / (r21^pf GCI21) differs from 1 by more than 10 %)")

    return GCIResult(
        phi=(phi1, phi2, phi3), r21=r21, r32=r32, p=p, s=s,
        phi_ext21=phi_ext21, phi_ext32=phi_ext32,
        e_a21=e_a21, e_a32=e_a32, e_ext21=e_ext21,
        gci_fine21=gci21, gci_medium32=gci32,
        asymptotic_ratio=asym, iterations=iters, warnings=warnings,
    )


def gci_from_cells(phi: Sequence[float], cells: Sequence[float], dim: int,
                   safety_factor: float = 1.25,
                   formal_order: Optional[float] = None) -> GCIResult:
    """Three-grid GCI from solution values and cell counts, finest first."""
    if len(phi) != 3 or len(cells) != 3:
        raise ValueError("give exactly three values and three cell counts, finest first")
    r21 = refinement_ratio(cells[0], cells[1], dim)
    r32 = refinement_ratio(cells[1], cells[2], dim)
    return gci_three_grid(phi[0], phi[1], phi[2], r21, r32, safety_factor, formal_order)


def gci_two_grid(phi1: float, phi2: float, r: float, p: float,
                 safety_factor: float = 3.0) -> float:
    """Roache's two-grid GCI on the fine grid with an assumed order p.

    Returns a fraction. Roache recommends a safety factor of 3 for two grids.
    """
    if r <= 1.0:
        raise ValueError("refinement ratio must be greater than 1")
    if phi1 == 0.0:
        raise ConvergenceError("relative error is undefined when phi1 is zero")
    return safety_factor * abs((phi1 - phi2) / phi1) / (r ** p - 1.0)


def profile_gci(phi1: Sequence[float], phi2: Sequence[float], phi3: Sequence[float],
                r21: float, r32: float, safety_factor: float = 1.25):
    """Point-by-point GCI along a profile sampled at the same locations on all grids.

    Points where the study is undefined (equal values, zero values) are returned
    as None. Also returns the average apparent order over the valid points, which
    Celik et al. suggest reporting as the global order.
    """
    if not (len(phi1) == len(phi2) == len(phi3)):
        raise ValueError("profiles must have the same length")
    results: List[Optional[GCIResult]] = []
    orders = []
    for a, b, c in zip(phi1, phi2, phi3):
        try:
            res = gci_three_grid(a, b, c, r21, r32, safety_factor)
        except ConvergenceError:
            results.append(None)
            continue
        results.append(res)
        orders.append(res.p)
    p_avg = sum(orders) / len(orders) if orders else float("nan")
    return results, p_avg
