import math

import pytest

from gcikit import (
    ConvergenceError,
    apparent_order,
    gci_from_cells,
    gci_three_grid,
    gci_two_grid,
    profile_gci,
    refinement_ratio,
)
from gcikit.cli import main


# Celik et al. (2008), Table 1, first column: N = 18000, 8000, 4500 (2D),
# phi = 6.063, 5.972, 5.863. Published: r21 = 1.5, r32 = 1.333, p = 1.53,
# phi_ext21 = 6.1685, e_a21 = 1.5 %, e_ext21 = 1.7 %, GCI_fine21 = 2.2 %.
CELIK_PHI = (6.063, 5.972, 5.863)
CELIK_CELLS = (18000, 8000, 4500)


def test_celik_table1_column1():
    res = gci_from_cells(CELIK_PHI, CELIK_CELLS, dim=2)
    assert res.r21 == pytest.approx(1.5)
    assert res.r32 == pytest.approx(4 / 3)
    assert res.p == pytest.approx(1.53, abs=0.01)
    assert res.phi_ext21 == pytest.approx(6.1685, abs=5e-4)
    assert 100 * res.e_a21 == pytest.approx(1.5, abs=0.05)
    assert 100 * res.e_ext21 == pytest.approx(1.7, abs=0.05)
    assert 100 * res.gci_fine21 == pytest.approx(2.2, abs=0.05)
    assert not res.oscillatory


def test_exact_second_order_data_recovers_p2():
    # phi(h) = 1 + 0.001 h^2 sampled at h = 1, 2, 4
    phi = [1 + 0.001 * h ** 2 for h in (1.0, 2.0, 4.0)]
    res = gci_three_grid(*phi, r21=2.0, r32=2.0, formal_order=2.0)
    assert res.p == pytest.approx(2.0, abs=1e-6)
    assert res.phi_ext21 == pytest.approx(1.0, abs=1e-12)
    assert res.asymptotic_ratio == pytest.approx(1.0, abs=0.01)
    assert res.warnings == []


def test_outside_asymptotic_range_warns():
    # observed order 1 but the scheme is formally second order
    phi = [1 + 0.01 * h for h in (1.0, 2.0, 4.0)]
    res = gci_three_grid(*phi, r21=2.0, r32=2.0, formal_order=2.0)
    assert res.p == pytest.approx(1.0, abs=1e-6)
    assert any("asymptotic" in w for w in res.warnings)
    assert gci_three_grid(*phi, r21=2.0, r32=2.0).asymptotic_ratio is None


def test_unequal_ratios_recover_order():
    # phi(h) = 2 - 0.3 h^1.7 with h = 1, 1.4, 2.1 (r21 = 1.4, r32 = 1.5)
    hs = (1.0, 1.4, 2.1)
    phi = [2 - 0.3 * h ** 1.7 for h in hs]
    res = gci_three_grid(*phi, r21=1.4, r32=1.5)
    assert res.p == pytest.approx(1.7, abs=1e-8)
    assert res.phi_ext21 == pytest.approx(2.0, abs=1e-10)


def test_oscillatory_convergence_flagged():
    res = gci_three_grid(1.00, 1.02, 0.97, r21=2.0, r32=2.0)
    assert res.s == -1
    assert res.oscillatory
    assert any("oscillatory" in w for w in res.warnings)


def test_identical_values_raise():
    with pytest.raises(ConvergenceError):
        gci_three_grid(1.0, 1.0, 1.1, r21=2.0, r32=2.0)


def test_bad_ratio_raises():
    with pytest.raises(ValueError):
        gci_three_grid(1.0, 1.1, 1.2, r21=1.0, r32=2.0)


def test_refinement_ratio_checks():
    assert refinement_ratio(8000, 1000, 3) == pytest.approx(2.0)
    with pytest.raises(ValueError):
        refinement_ratio(1000, 8000, 3)
    with pytest.raises(ValueError):
        refinement_ratio(8000, 1000, 4)


def test_apparent_order_equal_ratio_closed_form():
    # with r21 == r32 and s = 1, q = 0 and p = ln(eps32/eps21)/ln r
    p, s, _ = apparent_order(eps21=0.1, eps32=0.4, r21=2.0, r32=2.0)
    assert s == 1
    assert p == pytest.approx(math.log(4) / math.log(2))


def test_two_grid_roache():
    # fs = 3, r = 2, p = 2: GCI = 3 * |e| / 3 = |e|
    assert gci_two_grid(1.0, 1.01, r=2.0, p=2.0) == pytest.approx(0.01)


def test_profile_gci_skips_undefined_points():
    phi1 = [1 + 1.0, 5.0, 3 + 0.25]
    phi2 = [1 + 4.0, 5.0, 3 + 1.0]
    phi3 = [1 + 16.0, 5.5, 3 + 4.0]
    results, p_avg = profile_gci(phi1, phi2, phi3, r21=2.0, r32=2.0)
    assert results[1] is None
    assert results[0].p == pytest.approx(2.0)
    assert p_avg == pytest.approx(2.0)


def test_cli_text_and_json(capsys):
    assert main(["--phi", *map(str, CELIK_PHI), "--cells", *map(str, CELIK_CELLS),
                 "--dim", "2"]) == 0
    assert "GCI_fine21" in capsys.readouterr().out
    assert main(["--phi", "1", "1", "2", "--ratios", "2", "2"]) == 2
