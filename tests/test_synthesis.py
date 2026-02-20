"""Tests for audio synthesis."""

import numpy as np
import pytest
from rdkit import Chem
from chemsymphony.features import extract_all_features
from chemsymphony.config import load_config
from chemsymphony.mapping.master import apply_master_mapping
from chemsymphony.mapping import generate_all_layers
from chemsymphony.synthesis.audio_engine import synthesize_audio, synthesize_layer
from chemsymphony.synthesis.mixer import mix_layers
from chemsymphony.synthesis.effects import apply_reverb, apply_eq, fade_in_out


def _setup(smiles: str):
    mol = Chem.MolFromSmiles(smiles)
    feat = extract_all_features(mol)
    cfg = load_config()
    apply_master_mapping(feat, cfg)
    layers = generate_all_layers(feat, cfg)
    return feat, cfg, layers


class TestSynthesis:
    def test_synthesize_ethanol(self):
        feat, cfg, layers = _setup("CCO")
        raw = synthesize_audio(layers, feat, cfg)
        assert len(raw) > 0
        for name, audio, vol in raw:
            assert isinstance(audio, np.ndarray)
            assert len(audio) > 0

    def test_mix_produces_stereo(self):
        feat, cfg, layers = _setup("CCO")
        raw = synthesize_audio(layers, feat, cfg)
        mixed = mix_layers(raw, feat, cfg)
        assert mixed.ndim == 2
        assert mixed.shape[1] == 2

    def test_mix_normalized(self):
        feat, cfg, layers = _setup("c1ccccc1")
        raw = synthesize_audio(layers, feat, cfg)
        mixed = mix_layers(raw, feat, cfg)
        assert np.max(np.abs(mixed)) <= 1.0


class TestEffects:
    def test_reverb(self):
        sr = 44100
        audio = np.sin(2 * np.pi * 440 * np.arange(sr) / sr)
        result = apply_reverb(audio, sr, depth=0.3)
        assert len(result) == len(audio)

    def test_eq(self):
        sr = 44100
        audio = np.sin(2 * np.pi * 440 * np.arange(sr) / sr)
        result = apply_eq(audio, sr, brightness=0.5)
        assert len(result) == len(audio)

    def test_fade(self):
        sr = 44100
        audio = np.ones(sr)
        result = fade_in_out(audio, sr, fade_ms=10)
        assert result[0] < 0.5  # Faded in
        assert result[-1] < 0.5  # Faded out
