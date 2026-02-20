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
    """Generate percussion layers from halogen/element distribution.

    Uses swing for alternate-hit displacement, velocity accents on downbeats,
    and arrangement_density to gate percussion activity.
    """
    ap = feat.audio_parameters
    bpm = ap["bpm"]
    duration_sec = ap["duration"]
    beats_total = (duration_sec / 60.0) * bpm

    # New audio params
    swing = ap.get("swing", 0.0)
    arrangement_density = ap.get("arrangement_density", 1.0)

    layers: list[Layer] = []

    # Halogen percussion
    perc_idx = 0
    for element, perc_note in _ELEMENT_PERC_NOTE.items():
        count = feat.element_counts.get(element, 0)
        if count == 0:
            continue

        # Density gating: skip some percussion layers if density is low
        perc_idx += 1
        if arrangement_density < 0.8 and perc_idx > 1:
            continue  # Only keep the first percussion layer for sparse arrangements

        # Distribute hits evenly across the duration
        positions = feat.element_positions.get(element, [])
        n_atoms = max(len(positions), count)
        interval = beats_total / n_atoms if n_atoms > 0 else beats_total

        notes: list[NoteEvent] = []
        for k in range(n_atoms):
            cluster_size = feat.element_clustering.get(element, 1)
            base_beat = k * interval

            if cluster_size > 1:
                for sub in range(cluster_size):
                    hit_beat = base_beat + sub * 0.25
                    # Swing on alternate hits
                    if swing > 0 and sub % 2 == 1:
                        hit_beat += swing * 0.25

                    # Velocity accent on downbeats
                    vel = 100 if sub == 0 else 75
                    notes.append(NoteEvent(
                        pitch=perc_note,
                        start=hit_beat,
                        duration=0.25,
                        velocity=vel,
                        channel=9,
                    ))
            else:
                hit_beat = base_beat
                # Swing on alternate main hits
                if swing > 0 and k % 2 == 1:
                    hit_beat += swing * interval * 0.5

                # Velocity accent on downbeats (every 4th hit)
                vel = 100 if k % 4 == 0 else 80
                notes.append(NoteEvent(
                    pitch=perc_note,
                    start=hit_beat,
                    duration=0.25,
                    velocity=vel,
                    channel=9,
                ))

        layers.append(Layer(
            name=f"percussion_{element}",
            instrument=f"percussion_{element}",
            channel=9,
            notes=notes,
            program=0,
        ))

    # Universal heartbeat rhythm: all molecules with >5 heavy atoms get a
    # rhythmic foundation based on atomic properties.
    if feat.heavy_atom_count > 5:
        heartbeat_notes: list[NoteEvent] = []

        # Ring count determines time feel
        if feat.ring_count == 0:
            beats_per_measure = 4  # 4/4
        elif feat.ring_count <= 2:
            beats_per_measure = 3  # 3/4
        else:
            beats_per_measure = 6  # 6/8

        # Double bond fraction determines subdivision density
        subdivisions = 1
        if feat.double_bond_fraction > 0.3:
            subdivisions = 2
        if feat.double_bond_fraction > 0.6:
            subdivisions = 4

        # Aromatic fraction modulates velocity
        base_velocity = int(60 + feat.aromatic_fraction * 40)
        base_velocity = min(100, base_velocity)

        beat = 0.0
        step_size = 1.0 / subdivisions
        measure_beat = 0
        while beat < beats_total:
            # Downbeat accent
            is_downbeat = measure_beat == 0
            vel = min(120, base_velocity + 20) if is_downbeat else base_velocity

            # Apply swing to off-beat subdivisions
            hit_beat = beat
            if swing > 0 and measure_beat % 2 == 1:
                hit_beat += swing * step_size * 0.5

            heartbeat_notes.append(NoteEvent(
                pitch=36,  # Bass Drum 1
                start=hit_beat,
                duration=0.25,
                velocity=vel,
                channel=9,
            ))

            beat += step_size
            measure_beat = (measure_beat + 1) % (beats_per_measure * subdivisions)

        if heartbeat_notes:
            layers.append(Layer(
                name="percussion_heartbeat",
                instrument="percussion_heartbeat",
                channel=9,
                notes=heartbeat_notes,
                program=0,
            ))

    return layers
