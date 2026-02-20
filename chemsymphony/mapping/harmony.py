"""§6 mapping: Branches → counter-melodies & harmonies."""

from __future__ import annotations

from chemsymphony.config import Config
from chemsymphony.features import MolecularFeatures
from chemsymphony.mapping.melody import Layer, NoteEvent

# Element → GM program number
_ELEMENT_PROGRAMS: dict[str, int] = {
    "C": 0,    # Piano
    "O": 73,   # Flute
    "N": 60,   # French Horn
    "S": 38,   # Synth Bass
    "P": 12,   # Marimba
    "F": 0, "Cl": 0, "Br": 0, "I": 0,
}

# Depth → interval offsets from root
_DEPTH_INTERVALS = {
    1: [4, 7],      # Third and fifth
    2: [10, 14],    # Seventh and ninth
    3: [6, 11],     # Tritone and major seventh
}


def generate_harmony(feat: MolecularFeatures, cfg: Config) -> list[Layer]:
    """Generate counter-melody layers — one per branch."""
    if feat.branch_count == 0:
        return []

    ap = feat.audio_parameters
    root = ap["root_note_midi"]
    scale = ap["scale_intervals"]
    bpm = ap["bpm"]
    duration_sec = ap["duration"]
    beats_total = (duration_sec / 60.0) * bpm
    chain_len = max(feat.longest_chain, 1)

    layers: list[Layer] = []

    for i, branch in enumerate(feat.branches):
        length = branch["length"]
        depth = branch["depth"]
        chain_pos = branch["chain_position"]
        symbols = branch["symbols"]
        is_symmetric = branch.get("is_symmetric", False)

        # Entry time based on chain position
        entry_beat = (chain_pos / chain_len) * beats_total

        # Intervals based on depth
        intervals = _DEPTH_INTERVALS.get(min(depth, 3), [4, 7])

        # Determine instrument from most common element in branch
        if symbols:
            from collections import Counter
            most_common = Counter(symbols).most_common(1)[0][0]
            program = _ELEMENT_PROGRAMS.get(most_common, 0)
        else:
            program = 0

        # Generate notes
        notes: list[NoteEvent] = []
        note_dur = max(0.25, min(2.0, beats_total / max(length * 4, 1)))

        if length <= 2:
            # Grace note or trill
            for k in range(length):
                interval = intervals[k % len(intervals)]
                pitch = root + interval
                notes.append(NoteEvent(
                    pitch=max(36, min(96, pitch)),
                    start=entry_beat + k * 0.25,
                    duration=0.25,
                    velocity=70,
                    channel=min(2 + i, 15),
                ))
        else:
            # Full counter-melodic phrase
            scale_deg = 0
            direction = 1
            for k in range(length):
                octave = scale_deg // len(scale)
                deg = scale_deg % len(scale)
                base_interval = intervals[0] if intervals else 4
                pitch = root + scale[deg] + 12 * octave + base_interval
                notes.append(NoteEvent(
                    pitch=max(36, min(96, pitch)),
                    start=entry_beat + k * note_dur,
                    duration=note_dur * 0.9,
                    velocity=75,
                    channel=min(2 + i, 15),
                ))
                scale_deg += direction
                if scale_deg >= len(scale) * 2 or scale_deg < 0:
                    direction *= -1

        # Symmetric branches: unison (no pan offset)
        pan = 0.0 if is_symmetric else (-0.5 + i * 0.3)
        pan = max(-1.0, min(1.0, pan))
        for note in notes:
            note.pan = pan

        layers.append(Layer(
            name=f"counter_melody_{i + 1}",
            instrument="counter_melody",
            channel=min(2 + i, 15),
            notes=notes,
            program=program,
        ))

    return layers
