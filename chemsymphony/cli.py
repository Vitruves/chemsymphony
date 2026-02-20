"""CLI for ChemSymphony using argparse with rich-argparse colored help."""

from __future__ import annotations

import argparse
import sys

from rich_argparse import RichHelpFormatter

from chemsymphony import __version__
from chemsymphony.canonicalize import InvalidSmilesError
from chemsymphony.config import load_config
from chemsymphony.core import ChemSymphony


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="chemsymphony",
        description="Generate unique audio from molecular SMILES strings.",
        formatter_class=RichHelpFormatter,
    )

    # Required
    parser.add_argument(
        "-s", "--smiles", required=True,
        help="Input SMILES string",
    )

    # Output options
    output_group = parser.add_argument_group("output options")
    output_group.add_argument(
        "-o", "--output", default=None,
        help="Output file or directory (default: ./output.<format>)",
    )
    output_group.add_argument(
        "-f", "--format", dest="fmt", default="mp3",
        choices=["mp3", "wav", "midi"],
        help="Output format (default: mp3)",
    )

    # Audio tuning
    tuning_group = parser.add_argument_group("audio tuning")
    tuning_group.add_argument(
        "--bpm", type=int, default=None,
        help="Override base tempo",
    )
    tuning_group.add_argument(
        "--duration", type=float, default=None,
        help="Override duration in seconds",
    )
    tuning_group.add_argument(
        "--key", default=None,
        help='Force musical key, e.g. "Cm", "G"',
    )
    tuning_group.add_argument(
        "--seed", type=int, default=None,
        help="Random seed for reproducibility",
    )

    # Verbosity
    verb_group = parser.add_argument_group("verbosity")
    verb_group.add_argument(
        "-v", "--verbose", action="store_true",
        help="Print extracted features and mapping details",
    )
    verb_group.add_argument(
        "-q", "--quiet", action="store_true",
        help="Suppress all non-error output",
    )
    verb_group.add_argument(
        "--dry-run", action="store_true",
        help="Extract features and print mapping, no audio",
    )

    parser.add_argument(
        "--version", action="version",
        version=f"chemsymphony {__version__}",
    )

    return parser


def main(argv: list[str] | None = None) -> None:
    """Entry point for the CLI."""
    parser = _build_parser()
    args = parser.parse_args(argv)

    cfg = load_config(bpm=args.bpm, duration=args.duration, key=args.key, seed=args.seed)
    cs = ChemSymphony(config=cfg)

    try:
        if args.dry_run:
            features = cs.extract_features(args.smiles)
            print(features.pretty())
            return

        output = args.output if args.output is not None else f"output.{args.fmt}"

        if args.verbose and not args.quiet:
            features = cs.extract_features(args.smiles)
            print(features.pretty())
            print()

        cs.generate_to_file(args.smiles, output=output, fmt=args.fmt)

        if not args.quiet:
            print(f"Wrote {args.fmt.upper()} to {output}")

    except InvalidSmilesError as exc:
        print(f"Error: {exc}", file=sys.stderr)
        sys.exit(1)
