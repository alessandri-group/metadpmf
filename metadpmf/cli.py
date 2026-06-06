"""metadpmf command-line interface."""

import argparse
import sys


def main():
    parser = argparse.ArgumentParser(
        prog="metadpmf",
        description="Metadynamics PMF workflow for GROMACS + PLUMED",
    )
    parser.add_argument(
        "-c", "--config",
        default="config.yaml",
        help="Path to config.yaml (default: ./config.yaml)",
    )

    subparsers = parser.add_subparsers(dest="command", metavar="COMMAND")
    subparsers.required = True

    # run: write plumed.dat + run_metad.sh
    subparsers.add_parser(
        "run",
        help="Write plumed.dat and run_metad.sh (grompp + mdrun)",
    )

    # reweight: write analysis/plumed_analysis.dat + reweight.sh
    subparsers.add_parser(
        "reweight",
        help="Write plumed_analysis.dat and reweight.sh (plumed driver)",
    )

    # fes: block FES analysis
    p = subparsers.add_parser(
        "fes",
        help="Run block FES analysis on analysis/COLVAR",
    )
    p.add_argument(
        "--block-size",
        type=int,
        default=None,
        metavar="N",
        help="Use a single block size instead of scanning (writes fes.dat directly)",
    )
    p.add_argument(
        "--1d",
        dest="marginal",
        action="store_true",
        help="2D run only: marginalise the 2D-biased COLVAR onto distance and "
             "write the 1D PMF (fes_1d.dat)",
    )

    # pmf: Jacobian correction + plot
    p = subparsers.add_parser(
        "pmf",
        help="Apply Jacobian correction, shift PMF to zero, and plot",
    )
    p.add_argument(
        "--temp", type=float, default=None,
        help="Override temperature (K)",
    )
    p.add_argument(
        "--1d",
        dest="marginal",
        action="store_true",
        help="2D run only: plot the 1D PMF from fes_1d.dat (run 'fes --1d' first)",
    )

    args = parser.parse_args()

    from metadpmf import config as cfg_mod

    try:
        cfg = cfg_mod.load(args.config)
    except (FileNotFoundError, ValueError) as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

    if args.command == "run":
        from metadpmf.steps import run
        run.run(cfg)

    elif args.command == "reweight":
        from metadpmf.steps import reweight
        reweight.run(cfg)

    elif args.command == "fes":
        from metadpmf.steps import fes
        fes.run(cfg, block_size=args.block_size, marginal=args.marginal)

    elif args.command == "pmf":
        from metadpmf.steps import pmf
        temp = args.temp or cfg["temperature"]
        pmf.run(cfg, temperature=temp, marginal=args.marginal)


if __name__ == "__main__":
    main()
