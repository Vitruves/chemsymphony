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

    # Required (unless --list-genres)
    parser.add_argument(
        "-s", "--smiles", default=None,
        help="Input SMILES string",
    )
    parser.add_argument(
        "--list-genres", action="store_true",
        help="List available genre profiles for --post-processing and exit",
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
    tuning_group.add_argument(
        "--post-processing", default=None, metavar="GENRE",
        help="Apply genre-specific audio styling (e.g. techno, ambient, gabber)",
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

    # --list-genres: print available genres and exit
    if args.list_genres:
        from chemsymphony.genres import list_genres, load_genre
        genres = list_genres()
        if not genres:
            print("No genres found.")
        else:
            print("Available genres for --post-processing:\n")
            for name in genres:
                profile = load_genre(name)
                bpm = f"{profile.bpm_range[0]}-{profile.bpm_range[1]} BPM" if profile.bpm_range else "any BPM"
                print(f"  {name:20s} {profile.name:20s} {bpm}")
        return

    # --smiles is required for generation
    if args.smiles is None:
        parser.error("the following arguments are required: -s/--smiles")

    # Validate genre name if provided
    if args.post_processing is not None:
        from chemsymphony.genres import list_genres, load_genre
        available = list_genres()
        if args.post_processing not in available:
            print(
                f"Error: unknown genre {args.post_processing!r}. "
                f"Available: {', '.join(available)}",
                file=sys.stderr,
            )
            sys.exit(1)

    cfg = load_config(
        bpm=args.bpm, duration=args.duration, key=args.key,
        seed=args.seed, post_processing=args.post_processing,
    )
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
