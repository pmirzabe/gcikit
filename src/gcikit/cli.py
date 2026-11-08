"""Command line: gcikit --phi F M C --cells NF NM NC --dim 3"""
from __future__ import annotations

import argparse
import json
import sys

from .gci import ConvergenceError, gci_from_cells, gci_three_grid


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        prog="gcikit",
        description="Three-grid GCI (Celik et al. 2008). Give values finest grid first.")
    ap.add_argument("--phi", nargs=3, type=float, required=True,
                    metavar=("FINE", "MEDIUM", "COARSE"), help="solution values")
    grp = ap.add_mutually_exclusive_group(required=True)
    grp.add_argument("--cells", nargs=3, type=float, metavar=("NF", "NM", "NC"),
                     help="cell counts, finest first (needs --dim)")
    grp.add_argument("--ratios", nargs=2, type=float, metavar=("R21", "R32"),
                     help="refinement ratios h2/h1 and h3/h2")
    ap.add_argument("--dim", type=int, choices=(1, 2, 3), help="problem dimension for --cells")
    ap.add_argument("--fs", type=float, default=1.25, help="safety factor (default 1.25)")
    ap.add_argument("--formal-order", type=float, default=None,
                    help="formal order of the scheme, for the asymptotic-range check")
    ap.add_argument("--json", action="store_true", help="print the result as JSON")
    args = ap.parse_args(argv)

    try:
        if args.cells:
            if args.dim is None:
                ap.error("--cells needs --dim")
            res = gci_from_cells(args.phi, args.cells, args.dim, args.fs, args.formal_order)
        else:
            res = gci_three_grid(*args.phi, *args.ratios, safety_factor=args.fs,
                                 formal_order=args.formal_order)
    except (ConvergenceError, ValueError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2

    if args.json:
        out = dict(res.__dict__)
        out["oscillatory"] = res.oscillatory
        print(json.dumps(out, indent=2))
    else:
        print(res.summary())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
