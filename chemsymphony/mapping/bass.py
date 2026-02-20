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
        size = ring["size"]
        heteroatoms = ring["heteroatoms"]
        is_aromatic = ring["is_aromatic"]
        is_fused = idx in fused_set
        substituent_count = ring.get("substituent_count", 0)

        # Base pitch with octave offset for ring system
        # Aromatic rings → higher register, saturated → lower
        octave_offset = ring_system_octave.get(idx, 0) * 12
        if is_aromatic:
            base_pitch = max(30, min(54, root - 12 + octave_offset))
        else:
            base_pitch = max(24, min(48, root - 24 + octave_offset))

        # Ring-size-based pitch patterns
        if size == 3:
            pitches = [base_pitch, base_pitch + 4, base_pitch + 7]  # Major triad (tension)
        elif size == 4:
            pitches = [base_pitch, base_pitch + 5, base_pitch, base_pitch + 5]  # Oscillating 4ths
        elif size == 5:
            pitches = [base_pitch, base_pitch + 3]  # Minor 3rd
        elif size == 7:
            pitches = [base_pitch, base_pitch + 5]  # Quartal (perfect 4th)
        elif size >= 8:
            pitches = [base_pitch, base_pitch + 3, base_pitch + 7, base_pitch + 10]  # Minor 7th arpeggio
        else:
            pitches = [base_pitch, base_pitch + 7]  # Perfect 5th (default, incl 6-ring)

        # Add chromatic color from heteroatoms
        for ha in heteroatoms:
            offset = _HETEROATOM_OFFSETS.get(ha, 2)
            pitches.append(base_pitch + offset)

        # Note duration and velocity based on ring character and topology
        # Fused rings get longer legato notes, isolated get staccato
        if is_fused:
            note_dur = 1.5 if is_aromatic else 1.0
        else:
            note_dur = 0.75 if is_aromatic else 0.35

        # Warmer filter = louder bass
        velocity_base = int(75 + filter_warmth * 30) if is_aromatic else int(90 + filter_warmth * 20)
        velocity_base = min(120, velocity_base)

        # Syncopation: heavily substituted rings (>=3) skip beats
        syncopated = substituent_count >= 3

        notes: list[NoteEvent] = []
        beat = 0.0
        step_count = min(size, 7)  # Cap loop length
        while beat < beats_total:
            for step in range(step_count):
                if beat + step >= beats_total:
                    break

                # Syncopated pattern: skip every other even step
                if syncopated and step % 3 == 1:
                    continue

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
                if step < step_count - 1 and not is_aromatic:
                    ghost_beat = start_beat + 0.5
                    if ghost_beat < beats_total:
                        notes.append(NoteEvent(
                            pitch=pitch,
                            start=ghost_beat,
                            duration=0.25,
                            velocity=max(30, velocity_base - 40),
                            channel=1,
                        ))

            beat += step_count

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
