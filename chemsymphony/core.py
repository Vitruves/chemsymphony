"""ChemSymphony main class — public API."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from chemsymphony.canonicalize import canonicalize
from chemsymphony.config import Config, load_config
from chemsymphony.features import MolecularFeatures, extract_all_features


@dataclass
class MidiResult:
    """Result of MIDI generation."""

    manifest: dict[str, Any]
    master: bytes
    tracks: list[tuple[str, bytes]]

    def save(self, directory: str | Path) -> None:
        """Save MIDI tracks + manifest to a directory."""
        from chemsymphony.export.midi_export import save_midi_directory
        save_midi_directory(self, Path(directory))


class ChemSymphony:
    """Main entry point for molecular audio generation."""

    def __init__(self, config: Config | None = None) -> None:
        self.config = config or load_config()

    def extract_features(self, smiles: str) -> MolecularFeatures:
        """Extract all molecular features without generating audio."""
        result = canonicalize(smiles)
        feat = extract_all_features(result.mol, result.canonical_smiles)
        feat.molecular_formula = feat.molecular_formula or ""

        from chemsymphony.mapping.master import apply_master_mapping
        apply_master_mapping(feat, self.config)

        return feat

    def generate(
        self,
        smiles: str,
        *,
        fmt: str = "mp3",
        bpm: int | None = None,
        key: str | None = None,
        seed: int | None = None,
    ) -> MidiResult | bytes:
        """Generate audio or MIDI from a SMILES string.

        Returns ``MidiResult`` for MIDI format, raw audio bytes otherwise.
        """
        cfg = self.config
        if bpm is not None:
            cfg.bpm = bpm
        if key is not None:
            cfg.key = key
        if seed is not None:
            cfg.seed = seed

        result = canonicalize(smiles)
        feat = extract_all_features(result.mol, result.canonical_smiles)

        from chemsymphony.mapping.master import apply_master_mapping
        apply_master_mapping(feat, cfg)

        from chemsymphony.mapping import generate_all_layers
        layers = generate_all_layers(feat, cfg)

        if fmt == "midi":
            from chemsymphony.synthesis.midi_builder import build_midi
            from chemsymphony.export.manifest import build_manifest
            midi_data = build_midi(layers, feat, cfg)
            manifest = build_manifest(result, feat, layers, cfg)
            return MidiResult(
                manifest=manifest,
                master=midi_data["master"],
                tracks=midi_data["tracks"],
            )
        else:
            from chemsymphony.synthesis.audio_engine import synthesize_audio
            from chemsymphony.synthesis.mixer import mix_layers
            from chemsymphony.export.audio_export import export_audio
            raw_layers = synthesize_audio(layers, feat, cfg)
            mixed = mix_layers(raw_layers, feat, cfg)
            return export_audio(mixed, fmt=fmt, config=cfg)

    def generate_to_file(
        self,
        smiles: str,
        *,
        output: str | Path = "output.mp3",
        fmt: str = "mp3",
        bpm: int | None = None,
        key: str | None = None,
        seed: int | None = None,
    ) -> None:
        """Generate audio/MIDI and save to *output*."""
        result = self.generate(smiles, fmt=fmt, bpm=bpm, key=key, seed=seed)
        output = Path(output)

        if isinstance(result, MidiResult):
            result.save(output)
        else:
            output.parent.mkdir(parents=True, exist_ok=True)
            output.write_bytes(result)
