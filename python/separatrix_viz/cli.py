from __future__ import annotations

import argparse
import json

from .io import load_table, prepare_field, separatrix_summary
from .plotting import plot_networks
from .compare import compare_files


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="sepviz", description="Visualize and inspect Separatrix Toolkit output")
    sub = parser.add_subparsers(dest="command", required=True)

    p_plot = sub.add_parser("plot", help="plot separatrix networks")
    p_plot.add_argument("separatrices")
    p_plot.add_argument("--saddles", default=None)
    p_plot.add_argument("--field", default=None)
    p_plot.add_argument("--layer", type=int, default=None)
    p_plot.add_argument("--output", default="separatrices.png")
    p_plot.add_argument("--formats", nargs="*", default=None, choices=["png", "pdf", "svg"])
    p_plot.add_argument("--dpi", type=int, default=300)
    p_plot.add_argument("--levels", type=int, default=25)
    p_plot.add_argument("--title", default=None)
    p_plot.add_argument("--no-annotations", action="store_true")
    p_plot.add_argument("--no-equal-aspect", action="store_true")
    p_plot.add_argument("--show", action="store_true")

    p_prepare = sub.add_parser("prepare", help="extract x, y, psi from a simulation table")
    p_prepare.add_argument("source")
    p_prepare.add_argument("output")
    p_prepare.add_argument("--x-col", type=int, default=1)
    p_prepare.add_argument("--y-col", type=int, default=2)
    p_prepare.add_argument("--psi-col", type=int, default=3)

    p_compare = sub.add_parser("compare", help="compare candidate and reference curves")
    p_compare.add_argument("reference")
    p_compare.add_argument("candidate")
    p_compare.add_argument("--reference-schema", choices=["xy", "toolkit"], default="xy")
    p_compare.add_argument("--candidate-schema", choices=["xy", "toolkit"], default="toolkit")
    p_compare.add_argument("--json", action="store_true")

    p_inspect = sub.add_parser("inspect", help="summarize a separatrices.dat file")
    p_inspect.add_argument("separatrices")
    p_inspect.add_argument("--json", action="store_true")

    return parser


def main() -> None:
    args = build_parser().parse_args()

    if args.command == "plot":
        outputs = plot_networks(
            args.separatrices,
            saddles=args.saddles,
            field=args.field,
            layer=args.layer,
            title=args.title,
            output=args.output,
            formats=args.formats,
            dpi=args.dpi,
            contour_levels=args.levels,
            annotate_saddles=not args.no_annotations,
            equal_aspect=not args.no_equal_aspect,
            show=args.show,
        )
        for path in outputs:
            print(f"Saved: {path}")
        return

    if args.command == "prepare":
        prepare_field(args.source, args.output, x_col=args.x_col, y_col=args.y_col, psi_col=args.psi_col)
        print(f"Saved: {args.output}")
        return

    if args.command == "compare":
        result = compare_files(
            args.reference,
            args.candidate,
            reference_schema=args.reference_schema,
            candidate_schema=args.candidate_schema,
        )
        if args.json:
            print(json.dumps(result, indent=2))
        else:
            for key, value in result.items():
                print(f"{key}: {value:.8e}")
        return

    if args.command == "inspect":
        summary = separatrix_summary(load_table(args.separatrices))
        if args.json:
            print(json.dumps(summary, indent=2))
        else:
            for key, value in summary.items():
                print(f"{key}: {value}")
        return


if __name__ == "__main__":
    main()
