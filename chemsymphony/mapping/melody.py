"""§5 mapping: Main carbon chain → lead melody."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING

from chemsymphony.config import Config
from chemsymphony.features import MolecularFeatures

if TYPE_CHECKING:
    from chemsymphony.mapping.form import FormParams


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


def generate_melody(feat: MolecularFeatures, cfg: Config,
                    form: FormParams | None = None) -> Layer:
    """Generate the lead melody layer from the longest chain.

    Uses harmonic_tension for chromatic passing tones, swing for rhythmic
    displacement, and a velocity contour (crescendo → diminuendo).
    Optionally uses FormParams for climax position and octave range.
    """
    ap = feat.audio_parameters
    root = ap["root_note_midi"]
    scale = ap["scale_intervals"]
    bpm = ap["bpm"]
    duration_sec = ap["duration"]
    chain_len = feat.longest_chain

    if chain_len == 0:
        return Layer(name="lead_melody", instrument="acoustic_piano", channel=0)

    # New audio params
    harmonic_tension = ap.get("harmonic_tension", 0.0)
    swing = ap.get("swing", 0.0)

    # Note duration: spread chain notes across the composition
    beats_total = (duration_sec / 60.0) * bpm
    note_dur = max(0.25, beats_total / chain_len)

    notes: list[NoteEvent] = []
    direction = 1  # 1 = ascending, -1 = descending
    scale_degree = 0
    branch_set = set(feat.branch_points_on_chain)
    hetero_set = set(feat.chain_heteroatom_positions)
    bond_orders = feat.chain_bond_orders

    # Climax position for velocity contour (from form params or default 0.6)
    climax_pos = form.climax_position if form is not None else 0.6

    # Octave range from form params — clamp pitch range accordingly
    octave_range = form.octave_range if form is not None else 2
    pitch_lo = max(36, root - 12)
    pitch_hi = min(96, root + 12 * octave_range)

    for i in range(chain_len):
        # At branch points, reverse direction
        chain_atom = feat.longest_chain_atoms[i] if i < len(feat.longest_chain_atoms) else i
        if chain_atom in branch_set:
            direction *= -1

        # Compute pitch
        is_chromatic = False
        if i in hetero_set:
            # Chromatic passing tone
            pitch = root + scale[scale_degree % len(scale)] + 1
            is_chromatic = True
        elif harmonic_tension > 0.3 and i % 4 == 3:
            # Occasional chromatic passing tone based on harmonic tension
            # Every 4th note gets chromaticism if tension is high enough
            degree_in_scale = scale_degree % len(scale)
            base_pitch = root + scale[degree_in_scale] + 12 * (scale_degree // len(scale))
            # Approach from a semitone below
            pitch = base_pitch - 1
            is_chromatic = True
        else:
            octave = scale_degree // len(scale)
            degree_in_scale = scale_degree % len(scale)
            pitch = root + scale[degree_in_scale] + 12 * octave

        pitch = max(pitch_lo, min(pitch_hi, pitch))

        # Determine note duration and effects based on bond orders
        dur = note_dur
        velocity = 90
        effects: dict = {}

        if i < len(bond_orders):
            bond_order = bond_orders[i]
            if bond_order == 2:
                dur = note_dur * 0.5
                velocity = 110
                effects["pitch_bend"] = True
            elif bond_order == 3:
                dur = note_dur * 2.0
                velocity = 85
                effects["vibrato"] = True

        # Velocity contour: crescendo toward climax, diminuendo after
        progress = i / max(chain_len - 1, 1)
        if progress <= climax_pos:
            # Crescendo: ramp from 0.7 to 1.0
            contour = 0.7 + 0.3 * (progress / climax_pos)
        else:
            # Diminuendo: ramp from 1.0 to 0.6
            contour = 1.0 - 0.4 * ((progress - climax_pos) / (1.0 - climax_pos))
        velocity = max(40, min(127, int(velocity * contour)))

        # Swing: displace off-beat notes
        start_beat = i * note_dur
        if swing > 0 and i % 2 == 1:
            # Off-beat notes get displaced forward
            start_beat += swing * note_dur

        notes.append(NoteEvent(
            pitch=pitch,
            start=start_beat,
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
