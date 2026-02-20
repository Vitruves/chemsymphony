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

        # Single sustained note spanning the composition
        notes = [NoteEvent(
            pitch=max(36, min(84, ring_root)),
            start=0.0,
            duration=beats_total * sustain_factor,
            velocity=pad_velocity,
            channel=4,
            effects=filter_effects,
        )]

        # Add a fifth above for richer texture
        notes.append(NoteEvent(
            pitch=max(36, min(84, ring_root + 7)),
            start=0.0,
            duration=beats_total * sustain_factor,
            velocity=pad_velocity - 10,
            channel=4,
            effects={"pad": True, "timbre_organic": timbre_organic, **detune_effect},
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
