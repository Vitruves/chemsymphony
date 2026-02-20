"""§2+11 mapping: Element composition → instrument palette & percussion patterns."""

from __future__ import annotations

from chemsymphony.config import Config
from chemsymphony.features import MolecularFeatures
from chemsymphony.mapping.melody import Layer, NoteEvent

# Element → GM percussion note (channel 9)
_ELEMENT_PERC_NOTE: dict[str, int] = {
    "F": 42,   # Closed Hi-Hat
    "Cl": 38,  # Acoustic Snare
    "Br": 49,  # Crash Cymbal
    "I": 53,   # Ride Bell (gong-like)
}

# Non-percussion elements → pitched GM program
_ELEMENT_PROGRAM: dict[str, int] = {
    "C": 0,     # Acoustic Grand Piano
    "H": 8,     # Celesta
    "O": 73,    # Flute
    "N": 60,    # French Horn
    "S": 38,    # Synth Bass 1
    "P": 12,    # Marimba
    "Na": 9,    # Glockenspiel
    "K": 9,     # Glockenspiel
    "Fe": 9,    # Glockenspiel
}


def generate_percussion(feat: MolecularFeatures, cfg: Config) -> list[Layer]:
    """Generate percussion layers from halogen/element distribution."""
    ap = feat.audio_parameters
    bpm = ap["bpm"]
    duration_sec = ap["duration"]
    beats_total = (duration_sec / 60.0) * bpm

    layers: list[Layer] = []

    # Halogen percussion
    for element, perc_note in _ELEMENT_PERC_NOTE.items():
        count = feat.element_counts.get(element, 0)
        if count == 0:
            continue

        # Distribute hits evenly across the duration
        positions = feat.element_positions.get(element, [])
        n_atoms = max(len(positions), count)
        interval = beats_total / n_atoms if n_atoms > 0 else beats_total

        notes: list[NoteEvent] = []
        for k in range(n_atoms):
            # Check clustering for grouped notes (tuplets)
            cluster_size = feat.element_clustering.get(element, 1)
            if cluster_size > 1:
                for sub in range(cluster_size):
                    notes.append(NoteEvent(
                        pitch=perc_note,
                        start=k * interval + sub * 0.25,
                        duration=0.25,
                        velocity=90,
                        channel=9,
                    ))
            else:
                notes.append(NoteEvent(
                    pitch=perc_note,
                    start=k * interval,
                    duration=0.25,
                    velocity=90,
                    channel=9,
                ))

        layers.append(Layer(
            name=f"percussion_{element}",
            instrument=f"percussion_{element}",
            channel=9,
            notes=notes,
            program=0,
        ))

    return layers
