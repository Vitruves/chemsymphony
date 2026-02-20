"""§4 mapping: Aromaticity → drone pads & reverb."""

from __future__ import annotations

from chemsymphony.config import Config
from chemsymphony.features import MolecularFeatures
from chemsymphony.mapping.melody import Layer, NoteEvent


def _voicing_for_heteroatoms(heteroatoms: list[str]) -> list[int]:
    """Choose voicing intervals based on heteroatom identity."""
    if not heteroatoms:
        return [0, 7]  # Power chord

    hetero_set = set(heteroatoms)
    n_hetero = len(heteroatoms)

    if n_hetero >= 3:
        return [0, 4, 8]          # Augmented
    if n_hetero == 2:
        return [0, 2, 7, 10]     # Sus2 + 7th
    # Single heteroatom — choose by identity
    if "N" in hetero_set:
        return [0, 4, 7]         # Major triad
    if "O" in hetero_set:
        return [0, 5, 7]         # Sus4
    if "S" in hetero_set:
        return [0, 3, 7]         # Minor triad
    return [0, 4, 7]             # Default: major triad


def generate_pads(feat: MolecularFeatures, cfg: Config) -> list[Layer]:
    """Generate drone pad layers from rings.

    Aromatic rings produce supersaw pads; saturated rings produce quieter
    sine drone pads. Voicings depend on heteroatom identity.
    """
    if feat.ring_count == 0:
        return []

    ap = feat.audio_parameters
    root = ap["root_note_midi"]
    bpm = ap["bpm"]
    duration_sec = ap["duration"]
    beats_total = (duration_sec / 60.0) * bpm

    timbre_organic = ap.get("timbre_organic", 0.5)
    filter_warmth = ap.get("filter_warmth", 0.5)

    # Sustain from conjugation length
    sustain_factor = min(2.0, feat.conjugation_length / 6.0) if feat.conjugation_length > 0 else 0.5

    layers: list[Layer] = []
    pad_idx = 0

    for ring in feat.rings:
        is_aromatic = ring["is_aromatic"]
        size = ring["size"]
        heteroatoms = ring["heteroatoms"]

        if is_aromatic:
            ring_root = root - 12  # One octave below melody
            pad_velocity = int(40 + feat.aromatic_fraction * 60)
            pad_velocity = min(100, pad_velocity)
            pad_instrument = "synth_pad"
            pad_program = 89
            local_sustain = sustain_factor
            local_organic = timbre_organic
        else:
            # Saturated ring: quieter sine drone in lower register
            ring_root = root - 24
            pad_velocity = int(30 + feat.fsp3 * 30)
            pad_velocity = min(70, pad_velocity)
            pad_instrument = "synth_pad"
            pad_program = 89
            local_sustain = max(0.3, sustain_factor * 0.6)
            local_organic = 1.0  # Force organic/sine character

        hetero_count = len(heteroatoms)
        detune_effect = {"detune_cents": hetero_count * 8} if hetero_count > 0 else {}

        # Timbral richness widens the supersaw voice spread
        timbral_richness = ap.get("timbral_richness", 0.5)
        voice_spread = 0.5 + 0.5 * timbral_richness  # 0.5–1.0

        filter_effects = {
            "pad": True,
            "reverb": feat.aromatic_atom_count / 20.0 if is_aromatic else 0.1,
            "filter_cutoff_scale": 0.3 + 0.7 * (1.0 - filter_warmth),
            "timbre_organic": local_organic,
            "voice_spread": voice_spread,
            **detune_effect,
        }

        # Voicings based on heteroatom identity
        voicing_intervals = _voicing_for_heteroatoms(heteroatoms)

        # Stagger entry time
        entry_beat = pad_idx * (beats_total * 0.1)
        entry_beat = min(entry_beat, beats_total * 0.3)

        notes: list[NoteEvent] = []
        for vi, interval in enumerate(voicing_intervals):
            pitch = max(36, min(84, ring_root + interval))
            vel = max(30, pad_velocity - vi * 5)
            notes.append(NoteEvent(
                pitch=pitch,
                start=entry_beat,
                duration=beats_total * local_sustain,
                velocity=vel,
                channel=4,
                effects=filter_effects if vi == 0 else {
                    "pad": True, "timbre_organic": local_organic, **detune_effect,
                },
            ))

        layers.append(Layer(
            name=f"drone_pad_{pad_idx + 1}",
            instrument=pad_instrument,
            channel=4,
            notes=notes,
            program=pad_program,
        ))
        pad_idx += 1

    return layers
