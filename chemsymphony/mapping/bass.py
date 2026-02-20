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
    """Generate bass loop layers — one per ring."""
    if feat.ring_count == 0:
        return []

    ap = feat.audio_parameters
    root = ap["root_note_midi"] - 24  # Bass register, 2 octaves below
    bpm = ap["bpm"]
    duration_sec = ap["duration"]
    beats_total = (duration_sec / 60.0) * bpm

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
        sub_count = ring["substituent_count"]

        # Base pitch with octave offset for ring system
        octave_offset = ring_system_octave.get(idx, 0) * 12
        base_pitch = max(24, min(60, root + octave_offset))

        # Build loop pitches: root + fifth, plus chromatic from heteroatoms
        pitches = [base_pitch, base_pitch + 7]
        for ha in heteroatoms:
            offset = _HETEROATOM_OFFSETS.get(ha, 2)
            pitches.append(base_pitch + offset)

        # Generate loop notes, repeating to fill duration
        loop_beats = size
        note_dur = 1.0 if is_aromatic else 0.5  # Legato vs staccato
        velocity_base = 80 if is_aromatic else 95

        notes: list[NoteEvent] = []
        beat = 0.0
        while beat < beats_total:
            for step in range(loop_beats):
                if beat + step >= beats_total:
                    break
                pitch = pitches[step % len(pitches)]
                notes.append(NoteEvent(
                    pitch=pitch,
                    start=beat + step,
                    duration=note_dur,
                    velocity=velocity_base,
                    channel=1,
                ))
            beat += loop_beats

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
