"""§4 mapping: Aromaticity → drone pads & reverb."""

from __future__ import annotations

from chemsymphony.config import Config
from chemsymphony.features import MolecularFeatures
from chemsymphony.mapping.melody import Layer, NoteEvent


def generate_pads(feat: MolecularFeatures, cfg: Config) -> list[Layer]:
    """Generate drone pad layers — one per aromatic ring.

    Uses timbre_organic to choose between sine (organic) and saw (synthetic) pad,
    and filter_warmth to control pad filter cutoff in effects dict.
    """
    if feat.aromatic_ring_count == 0:
        return []

    ap = feat.audio_parameters
    root = ap["root_note_midi"]
    bpm = ap["bpm"]
    duration_sec = ap["duration"]
    beats_total = (duration_sec / 60.0) * bpm

    # New audio params
    timbre_organic = ap.get("timbre_organic", 0.5)
    filter_warmth = ap.get("filter_warmth", 0.5)

    # Volume from aromatic fraction
    pad_velocity = int(40 + feat.aromatic_fraction * 60)
    pad_velocity = min(100, pad_velocity)

    # Sustain from conjugation length
    sustain_factor = min(2.0, feat.conjugation_length / 6.0) if feat.conjugation_length > 0 else 0.5

    # Choose instrument based on timbre_organic
    # High organic (sp3) = sine pad, low organic = saw pad (supersaw)
    pad_instrument = "synth_pad"  # Both use _synth_pad which now does supersaw
    # But we pass filter info through effects
    pad_program = 89  # Pad 2 (warm)

    scale_intervals = ap.get("scale_intervals", [0, 2, 4, 5, 7, 9, 11])

    layers: list[Layer] = []
    aromatic_idx = 0

    for ring in feat.rings:
        if not ring["is_aromatic"]:
            continue

        size = ring["size"]
        ring_root = root - 12  # One octave below melody

        # Detuning from heteroatoms
        hetero_count = len(ring["heteroatoms"])
        detune_effect = {"detune_cents": hetero_count * 8} if hetero_count > 0 else {}

        # Pad filter cutoff driven by filter_warmth
        filter_effects = {
            "pad": True,
            "reverb": feat.aromatic_atom_count / 20.0,
            "filter_cutoff_scale": 0.3 + 0.7 * (1.0 - filter_warmth),
            "timbre_organic": timbre_organic,
            **detune_effect,
        }

        # Scale-aware voicings based on ring heteroatom count
        # 0 heteroatoms → power chord (root + 5th)
        # 1 → triad (root + 3rd + 5th)
        # 2 → seventh chord (root + 3rd + 5th + 7th)
        # 3+ → cluster chord (root + 2nd + 3rd + 5th)
        voicing_intervals = [0]  # Always include root
        if hetero_count == 0:
            # Power chord: root + fifth
            voicing_intervals.append(7)
        elif hetero_count == 1:
            # Triad from scale
            if len(scale_intervals) >= 5:
                voicing_intervals.append(scale_intervals[2])  # 3rd
                voicing_intervals.append(scale_intervals[4])  # 5th
            else:
                voicing_intervals.extend([4, 7])
        elif hetero_count == 2:
            # Seventh chord from scale
            if len(scale_intervals) >= 7:
                voicing_intervals.append(scale_intervals[2])  # 3rd
                voicing_intervals.append(scale_intervals[4])  # 5th
                voicing_intervals.append(scale_intervals[6])  # 7th
            else:
                voicing_intervals.extend([4, 7, 10])
        else:
            # Cluster chord
            if len(scale_intervals) >= 5:
                voicing_intervals.append(scale_intervals[1])  # 2nd
                voicing_intervals.append(scale_intervals[2])  # 3rd
                voicing_intervals.append(scale_intervals[4])  # 5th
            else:
                voicing_intervals.extend([2, 4, 7])

        # Stagger entry time for multiple aromatic rings
        entry_beat = aromatic_idx * (beats_total * 0.1)
        entry_beat = min(entry_beat, beats_total * 0.3)  # Don't start too late

        notes: list[NoteEvent] = []
        for vi, interval in enumerate(voicing_intervals):
            pitch = max(36, min(84, ring_root + interval))
            vel = max(40, pad_velocity - vi * 5)
            notes.append(NoteEvent(
                pitch=pitch,
                start=entry_beat,
                duration=beats_total * sustain_factor,
                velocity=vel,
                channel=4,
                effects=filter_effects if vi == 0 else {
                    "pad": True, "timbre_organic": timbre_organic, **detune_effect,
                },
            ))

        layers.append(Layer(
            name=f"drone_pad_{aromatic_idx + 1}",
            instrument=pad_instrument,
            channel=4,
            notes=notes,
            program=pad_program,
        ))
        aromatic_idx += 1

    return layers
