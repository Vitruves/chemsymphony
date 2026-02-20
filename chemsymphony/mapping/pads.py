"""§4 mapping: Aromaticity → drone pads & reverb."""

from __future__ import annotations

from chemsymphony.config import Config
from chemsymphony.features import MolecularFeatures
from chemsymphony.mapping.melody import Layer, NoteEvent


def generate_pads(feat: MolecularFeatures, cfg: Config) -> list[Layer]:
    """Generate drone pad layers — one per aromatic ring."""
    if feat.aromatic_ring_count == 0:
        return []

    ap = feat.audio_parameters
    root = ap["root_note_midi"]
    bpm = ap["bpm"]
    duration_sec = ap["duration"]
    beats_total = (duration_sec / 60.0) * bpm

    # Volume from aromatic fraction
    pad_velocity = int(40 + feat.aromatic_fraction * 60)
    pad_velocity = min(100, pad_velocity)

    # Sustain from conjugation length
    sustain_factor = min(2.0, feat.conjugation_length / 6.0) if feat.conjugation_length > 0 else 0.5

    layers: list[Layer] = []
    aromatic_idx = 0

    for ring in feat.rings:
        if not ring["is_aromatic"]:
            continue

        # Pad tuned to the ring's bass root
        size = ring["size"]
        ring_root = root - 12  # One octave below melody

        # Detuning from heteroatoms (±5-15 cents encoded as slight pitch shift)
        hetero_count = len(ring["heteroatoms"])
        detune_effect = {"detune_cents": hetero_count * 8} if hetero_count > 0 else {}

        # Single sustained note spanning the composition
        notes = [NoteEvent(
            pitch=max(36, min(84, ring_root)),
            start=0.0,
            duration=beats_total * sustain_factor,
            velocity=pad_velocity,
            channel=4,
            effects={"pad": True, "reverb": feat.aromatic_atom_count / 20.0, **detune_effect},
        )]

        # Add a fifth above for richer texture
        notes.append(NoteEvent(
            pitch=max(36, min(84, ring_root + 7)),
            start=0.0,
            duration=beats_total * sustain_factor,
            velocity=pad_velocity - 10,
            channel=4,
            effects={"pad": True, **detune_effect},
        ))

        layers.append(Layer(
            name=f"drone_pad_{aromatic_idx + 1}",
            instrument="synth_pad",
            channel=4,
            notes=notes,
            program=89,  # Pad 2 (warm)
        ))
        aromatic_idx += 1

    return layers
