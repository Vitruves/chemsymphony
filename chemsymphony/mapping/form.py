"""§10 mapping: Graph topology → musical structure & form."""

from __future__ import annotations

from dataclasses import dataclass, field

from chemsymphony.config import Config
from chemsymphony.features import MolecularFeatures


@dataclass
class FormParams:
    """High-level musical form parameters."""

    octave_range: int = 1
    polyphony: int = 1
    climax_position: float = 0.5  # 0–1 fraction of total duration
    section_count: int = 1
    is_palindrome: bool = False
    pause_positions: list[float] = field(default_factory=list)  # In beats


def generate_form(feat: MolecularFeatures, cfg: Config) -> FormParams:
    """Compute musical form parameters from graph topology."""
    form = FormParams()

    # Diameter → octave range
    d = feat.graph_diameter
    if d <= 3:
        form.octave_range = 1
    elif d <= 8:
        form.octave_range = 2
    else:
        form.octave_range = 3

    # Average degree → polyphony
    if feat.avg_degree < 2.0:
        form.polyphony = min(2, max(1, int(feat.avg_degree)))
    else:
        form.polyphony = min(6, int(feat.avg_degree * 1.5))

    # Max degree atom → climax position
    chain = feat.longest_chain_atoms
    if chain and feat.max_degree_atom in chain:
        pos = chain.index(feat.max_degree_atom)
        form.climax_position = pos / max(len(chain) - 1, 1)
    else:
        form.climax_position = 0.5

    # Connected components → sections
    form.section_count = feat.connected_components

    # Symmetry → palindrome
    form.is_palindrome = feat.symmetry_score > 0.3

    # Bridges → pause positions (in beats)
    ap = feat.audio_parameters
    bpm = ap.get("bpm", 120)
    duration_sec = ap.get("duration", 10.0)
    beats_total = (duration_sec / 60.0) * bpm
    chain_len = max(feat.longest_chain, 1)

    for a, b in feat.bridges:
        # Find if bridge atoms are on the chain
        if a in chain and b in chain:
            pos = chain.index(a)
            beat = (pos / chain_len) * beats_total
            form.pause_positions.append(beat)

    return form
