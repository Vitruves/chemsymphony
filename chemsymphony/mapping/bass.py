"""§3 mapping: Ring systems → bass loops & rhythmic layers."""

from __future__ import annotations

from chemsymphony.config import Config
from chemsymphony.features import MolecularFeatures
from chemsymphony.mapping.melody import Layer, NoteEvent

# Heteroatom → chromatic offset for bass loop pitches
_HETEROATOM_OFFSETS: dict[str, int] = {
    "N": 3,
    "O": 7,
    "S": 10,
    "P": 5,
    "F": 1,
    "Cl": 2,
    "Br": 4,
    "I": 6,
}


def generate_bass(feat: MolecularFeatures, cfg: Config) -> list[Layer]:
    """Generate bass loop layers — one per ring.

    Uses swing for off-beat displacement, filter_warmth for velocity,
    and adds ghost notes for groove.
    """
    if feat.ring_count == 0:
        return []

    ap = feat.audio_parameters
    root = ap["root_note_midi"] - 24  # Bass register, 2 octaves below
    bpm = ap["bpm"]
    duration_sec = ap["duration"]
    beats_total = (duration_sec / 60.0) * bpm

    # New audio params
    swing = ap.get("swing", 0.0)
    filter_warmth = ap.get("filter_warmth", 0.5)

    fused_set = set()
    for i, j in feat.fused_ring_pairs:
        fused_set.add(i)
        fused_set.add(j)

    layers: list[Layer] = []
    ring_system_octave: dict[int, int] = {}
    for sys_idx, system in enumerate(feat.ring_systems):
        for ring_idx in system:
            ring_system_octave[ring_idx] = sys_idx

    for ring in feat.rings:
        idx = ring["index"]
        size = min(ring["size"], 7)
        heteroatoms = ring["heteroatoms"]
        is_aromatic = ring["is_aromatic"]

        # Base pitch with octave offset for ring system
        octave_offset = ring_system_octave.get(idx, 0) * 12
        base_pitch = max(24, min(60, root + octave_offset))

        # Build loop pitches: root + fifth, plus chromatic from heteroatoms
        pitches = [base_pitch, base_pitch + 7]
        for ha in heteroatoms:
            offset = _HETEROATOM_OFFSETS.get(ha, 2)
            pitches.append(base_pitch + offset)

        # Note duration and velocity based on ring character
        note_dur = 1.0 if is_aromatic else 0.5
        # Warmer filter = louder bass
        velocity_base = int(75 + filter_warmth * 30) if is_aromatic else int(90 + filter_warmth * 20)
        velocity_base = min(120, velocity_base)

        notes: list[NoteEvent] = []
        beat = 0.0
        while beat < beats_total:
            for step in range(size):
                if beat + step >= beats_total:
                    break
                pitch = pitches[step % len(pitches)]

                # Apply swing to off-beat notes
                start_beat = beat + step
                if swing > 0 and step % 2 == 1:
                    start_beat += swing * note_dur

                notes.append(NoteEvent(
                    pitch=pitch,
                    start=start_beat,
                    duration=note_dur,
                    velocity=velocity_base,
                    channel=1,
                ))

                # Ghost note between main hits for groove
                if step < size - 1 and not is_aromatic:
                    ghost_beat = start_beat + 0.5
                    if ghost_beat < beats_total:
                        notes.append(NoteEvent(
                            pitch=pitch,
                            start=ghost_beat,
                            duration=0.25,
                            velocity=max(30, velocity_base - 40),
                            channel=1,
                        ))

            beat += size

        # Pan based on ring order
        total_rings = feat.ring_count
        pan = -0.8 + (1.6 * idx / max(total_rings - 1, 1)) if total_rings > 1 else 0.0
        for note in notes:
            note.pan = pan

        layers.append(Layer(
            name=f"bass_loop_ring{idx + 1}",
            instrument="synth_bass",
            channel=1,
            notes=notes,
            program=38,  # Synth Bass 1
        ))

    return layers
