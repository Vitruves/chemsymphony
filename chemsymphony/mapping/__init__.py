"""Mapping orchestrator — generates all audio layers from molecular features."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

from chemsymphony.config import Config
from chemsymphony.features import MolecularFeatures
from chemsymphony.mapping.melody import Layer, NoteEvent
from chemsymphony.mapping.expression import ExpressionParams
from chemsymphony.mapping.form import FormParams


@dataclass
class CompositionLayers:
    """All generated musical layers plus global parameters."""

    melody: Layer | None = None
    bass: list[Layer] = field(default_factory=list)
    harmony: list[Layer] = field(default_factory=list)
    pads: list[Layer] = field(default_factory=list)
    percussion: list[Layer] = field(default_factory=list)
    motifs: Layer | None = None
    expression: ExpressionParams | None = None
    form: FormParams | None = None

    def all_layers(self) -> list[Layer]:
        """Return a flat list of all Layer objects."""
        result: list[Layer] = []
        if self.melody:
            result.append(self.melody)
        result.extend(self.bass)
        result.extend(self.harmony)
        result.extend(self.pads)
        result.extend(self.percussion)
        if self.motifs:
            result.append(self.motifs)
        return result


def generate_all_layers(feat: MolecularFeatures, cfg: Config) -> CompositionLayers:
    """Run all mapping pipelines and return the complete layer set."""
    from chemsymphony.mapping.melody import generate_melody
    from chemsymphony.mapping.bass import generate_bass
    from chemsymphony.mapping.harmony import generate_harmony
    from chemsymphony.mapping.pads import generate_pads
    from chemsymphony.mapping.percussion import generate_percussion
    from chemsymphony.mapping.motifs import generate_motifs
    from chemsymphony.mapping.expression import generate_expression
    from chemsymphony.mapping.form import generate_form

    layers = CompositionLayers()
    # Compute form first so melody can use climax_position and octave_range
    layers.form = generate_form(feat, cfg)
    layers.melody = generate_melody(feat, cfg, form=layers.form)
    layers.bass = generate_bass(feat, cfg)
    layers.harmony = generate_harmony(feat, cfg)
    layers.pads = generate_pads(feat, cfg)
    layers.percussion = generate_percussion(feat, cfg)
    layers.motifs = generate_motifs(feat, cfg)
    layers.expression = generate_expression(feat, cfg)

    # Apply form-level adjustments
    _apply_form(layers, feat, cfg)
    _apply_expression(layers, feat, cfg)

    return layers


def _apply_form(layers: CompositionLayers, feat: MolecularFeatures, cfg: Config) -> None:
    """Apply form parameters (pauses, palindrome) to the melody."""
    if layers.form is None or layers.melody is None:
        return

    form = layers.form

    # Insert pauses at bridge positions
    for pause_beat in form.pause_positions:
        for note in layers.melody.notes:
            if abs(note.start - pause_beat) < 0.5:
                note.velocity = max(0, note.velocity - 40)
                note.duration *= 0.3

    # Palindrome: mirror the melody
    if form.is_palindrome and layers.melody.notes:
        original = list(layers.melody.notes)
        total_dur = max(n.start + n.duration for n in original) if original else 0
        mirrored = []
        for note in reversed(original):
            mirrored.append(NoteEvent(
                pitch=note.pitch,
                start=total_dur + (total_dur - note.start),
                duration=note.duration,
                velocity=note.velocity,
                channel=note.channel,
                pan=note.pan,
                effects=dict(note.effects),
            ))
        layers.melody.notes.extend(mirrored)

    # Section form application
    if layers.melody.notes:
        original = list(layers.melody.notes)
        total_dur = max(n.start + n.duration for n in original) if original else 0
        section_form = getattr(form, "section_form", "through")

        if section_form == "aba" and total_dur > 0:
            # Repeat first half transposed down an octave
            half_dur = total_dur / 2
            reprise = []
            for note in original:
                if note.start < half_dur:
                    reprise.append(NoteEvent(
                        pitch=max(36, note.pitch - 12),
                        start=total_dur + note.start,
                        duration=note.duration,
                        velocity=max(40, note.velocity - 10),
                        channel=note.channel,
                        pan=note.pan,
                        effects=dict(note.effects),
                    ))
            layers.melody.notes.extend(reprise)
        elif section_form == "abab" and total_dur > 0:
            # Repeat the entire melody
            repeat = []
            for note in original:
                repeat.append(NoteEvent(
                    pitch=note.pitch,
                    start=total_dur + note.start,
                    duration=note.duration,
                    velocity=max(40, note.velocity - 15),
                    channel=note.channel,
                    pan=note.pan,
                    effects=dict(note.effects),
                ))
            layers.melody.notes.extend(repeat)


def _apply_expression(layers: CompositionLayers, feat: MolecularFeatures, cfg: Config) -> None:
    """Apply stereo and expression modifiers to all layers."""
    if layers.expression is None:
        return

    expr = layers.expression

    # Apply stereo width to all notes
    for layer in layers.all_layers():
        for note in layer.notes:
            note.pan *= expr.stereo_width

    # Apply charge accents
    if layers.melody and layers.melody.notes:
        chain = feat.longest_chain_atoms
        chain_len = max(len(chain), 1)
        bpm = feat.audio_parameters.get("bpm", 120)
        duration_sec = feat.audio_parameters.get("duration", 10.0)
        beats_total = (duration_sec / 60.0) * bpm

        for accent in expr.charge_accents:
            atom_idx = accent["atom_idx"]
            if atom_idx in chain:
                pos = chain.index(atom_idx)
                target_beat = (pos / chain_len) * beats_total
                for note in layers.melody.notes:
                    if abs(note.start - target_beat) < 0.5:
                        if accent["accent"] == "sforzando":
                            note.velocity = min(127, note.velocity + 30)
                        else:
                            note.velocity = max(40, note.velocity - 20)
                            note.duration *= 1.5
