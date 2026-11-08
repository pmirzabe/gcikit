# gcikit

Discretization-uncertainty estimates for CFD and FEA mesh studies: the Grid
Convergence Index (GCI), apparent order of convergence and Richardson
extrapolation, following Celik et al. (2008) and Roache.

Most mesh-independence studies end in a spreadsheet. gcikit puts the standard
procedure in one tested Python package, so the numbers in a V&V report can be
reproduced from the solver output.

## Status

Early stage (v0.1). Working now:

- three-grid GCI with the apparent order solved by fixed-point iteration
  (Celik et al. 2008, eqs. 1-7), including oscillatory convergence;
- refinement ratios from cell counts in 1, 2 or 3 dimensions;
- Roache's two-grid GCI with an assumed order;
- point-by-point GCI along a profile, with the average apparent order;
- checks and warnings: refinement ratio below 1.3, unusual order, asymptotic range;
- a command-line tool.

The test suite reproduces the worked example in Table 1 of Celik et al. (2008).

## Install

```bash
pip install -e ".[dev]"
pytest
```

## Use

```python
from gcikit import gci_from_cells

res = gci_from_cells(phi=[6.063, 5.972, 5.863], cells=[18000, 8000, 4500], dim=2)
print(res.summary())
```

```bash
gcikit --phi 6.063 5.972 5.863 --cells 18000 8000 4500 --dim 2
```

Values are always given finest grid first.

## Roadmap

1. Readers for solver output: OpenFOAM `postProcessing/` function objects
   (forces, probes, surface field values) and residual logs, and ANSYS Fluent
   report files, so a study can be run straight from three case folders.
2. Least-squares GCI for more than three grids (Eça and Hoekstra 2014).
3. Profile and field studies: interpolation of profiles onto common points,
   and global and local orders.
4. Iterative-convergence checks, so iteration error is not mistaken for
   discretization error.
5. Plots (solution against h with the extrapolated value and error bands) and
   a Markdown or HTML report in the format Celik et al. recommend.
6. Worked examples on public cases: lid-driven cavity, backward-facing step.

## References

- I. B. Celik et al., "Procedure for Estimation and Reporting of Uncertainty Due
  to Discretization in CFD Applications", J. Fluids Eng. 130(7), 078001, 2008.
- P. J. Roache, "Verification and Validation in Computational Science and
  Engineering", Hermosa, 1998.
- L. Eça, M. Hoekstra, "A procedure for the estimation of the numerical
  uncertainty of CFD calculations based on grid refinement studies",
  J. Comput. Phys. 262, 2014.

## License

MIT
