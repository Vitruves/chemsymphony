"""§5 mapping: Main carbon chain → lead melody."""

from __future__ import annotations

from dataclasses import dataclass, field

from chemsymphony.config import Config
from chemsymphony.features import MolecularFeatures


@dataclass
class NoteEvent:
    """A single note event in any layer."""

    pitch: int          # MIDI note number (0–127)
    start: float        # Start time in beats
    duration: float     # Duration in beats
    velocity: int = 100 # MIDI velocity (0–127)
    channel: int = 0
    pan: float = 0.0    # -1.0 (left) to 1.0 (right)
    effects: dict = field(default_factory=dict)


@dataclass
class Layer:
    """A named collection of note events."""

    name: str
    instrument: str
    channel: int
    notes: list[NoteEvent] = field(default_factory=list)
    program: int = 0  # GM MIDI program number


def generate_melody(feat: MolecularFeatures, cfg: Config) -> Layer:
    """Generate the lead melody layer from the longest chain."""
    ap = feat.audio_parameters
    root = ap["root_note_midi"]
    scale = ap["scale_intervals"]
    bpm = ap["bpm"]
    duration_sec = ap["duration"]
    chain_len = feat.longest_chain

    if chain_len == 0:
        return Layer(name="lead_melody", instrument="acoustic_piano", channel=0)

    # Note duration: spread chain notes across the composition
    beats_total = (duration_sec / 60.0) * bpm
    note_dur = max(0.25, beats_total / chain_len)

    notes: list[NoteEvent] = []
    direction = 1  # 1 = ascending, -1 = descending
    scale_degree = 0
    branch_set = set(feat.branch_points_on_chain)
    hetero_set = set(feat.chain_heteroatom_positions)
    bond_orders = feat.chain_bond_orders

    for i in range(chain_len):
        # At branch points, reverse direction
        chain_atom = feat.longest_chain_atoms[i] if i < len(feat.longest_chain_atoms) else i
        if chain_atom in branch_set:
            direction *= -1

        # Compute pitch
        if i in hetero_set:
            # Chromatic passing tone
            pitch = root + scale[scale_degree % len(scale)] + 1
        else:
            octave = scale_degree // len(scale)
            degree_in_scale = scale_degree % len(scale)
            pitch = root + scale[degree_in_scale] + 12 * octave

        pitch = max(36, min(96, pitch))

        # Determine note duration and effects based on bond orders
        dur = note_dur
        velocity = 90
        effects: dict = {}

        if i < len(bond_orders):
            bond_order = bond_orders[i]
            if bond_order == 2:
                # Double bond: staccato + accent
                dur = note_dur * 0.5
                velocity = 110
                effects["pitch_bend"] = True
            elif bond_order == 3:
                # Triple bond: sustained + vibrato
                dur = note_dur * 2.0
                velocity = 85
                effects["vibrato"] = True

        notes.append(NoteEvent(
            pitch=pitch,
            start=i * note_dur,
            duration=dur,
            velocity=velocity,
            channel=0,
            effects=effects,
        ))

        scale_degree += direction

    return Layer(
        name="lead_melody",
        instrument="acoustic_piano",
        channel=0,
        notes=notes,
        program=0,
    )
