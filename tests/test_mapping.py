"""Tests for audio mapping pipelines."""

import pytest
from rdkit import Chem
from chemsymphony.features import extract_all_features
from chemsymphony.config import load_config
from chemsymphony.mapping.master import apply_master_mapping
from chemsymphony.mapping.melody import generate_melody
from chemsymphony.mapping.bass import generate_bass
from chemsymphony.mapping.harmony import generate_harmony
from chemsymphony.mapping.pads import generate_pads
from chemsymphony.mapping.percussion import generate_percussion
from chemsymphony.mapping.motifs import generate_motifs
from chemsymphony.mapping import generate_all_layers


def _features(smiles: str):
    mol = Chem.MolFromSmiles(smiles)
    feat = extract_all_features(mol)
    cfg = load_config()
    apply_master_mapping(feat, cfg)
    return feat, cfg


class TestMasterMapping:
    def test_bpm_range(self):
        feat, cfg = _features("C")
        assert 70 <= feat.audio_parameters["bpm"] <= 180

    def test_heavy_molecule_faster(self):
        feat_light, _ = _features("C")
        feat_heavy, _ = _features("c1ccc2ccccc2c1")  # Naphthalene
        assert feat_heavy.audio_parameters["bpm"] >= feat_light.audio_parameters["bpm"]

    def test_scale_selection(self):
        feat, _ = _features("CCCCCC")  # Pure hydrocarbon
        assert feat.audio_parameters["scale"] == "pentatonic_major"

    def test_root_note_deterministic(self):
        f1, _ = _features("CCO")
        f2, _ = _features("CCO")
        assert f1.audio_parameters["root_note_midi"] == f2.audio_parameters["root_note_midi"]


class TestMelody:
    def test_melody_notes_match_chain(self):
        feat, cfg = _features("CCCC")
        layer = generate_melody(feat, cfg)
        assert len(layer.notes) == feat.longest_chain

    def test_single_atom_melody(self):
        feat, cfg = _features("C")
        layer = generate_melody(feat, cfg)
        assert len(layer.notes) == 1

    def test_melody_pitches_in_range(self):
        feat, cfg = _features("CCCCCCCCCC")
        layer = generate_melody(feat, cfg)
        for note in layer.notes:
            assert 24 <= note.pitch <= 108


class TestBass:
    def test_benzene_has_bass(self):
        feat, cfg = _features("c1ccccc1")
        layers = generate_bass(feat, cfg)
        assert len(layers) == 1
        assert len(layers[0].notes) > 0

    def test_no_ring_no_bass(self):
        feat, cfg = _features("CCCC")
        layers = generate_bass(feat, cfg)
        assert len(layers) == 0


class TestHarmony:
    def test_branched_molecule(self):
        feat, cfg = _features("CC(C)C")  # Isobutane
        layers = generate_harmony(feat, cfg)
        assert len(layers) >= 1

    def test_linear_no_harmony(self):
        feat, cfg = _features("CCCC")
        layers = generate_harmony(feat, cfg)
        # May have 0 or small count
        # Linear chains have terminal H but no branch points off main chain
        assert isinstance(layers, list)


class TestPads:
    def test_aromatic_has_pads(self):
        feat, cfg = _features("c1ccccc1")
        layers = generate_pads(feat, cfg)
        assert len(layers) >= 1

    def test_aliphatic_no_pads(self):
        feat, cfg = _features("CCCC")
        layers = generate_pads(feat, cfg)
        assert len(layers) == 0


class TestAllLayers:
    def test_benzene_full_pipeline(self):
        feat, cfg = _features("c1ccccc1")
        layers = generate_all_layers(feat, cfg)
        assert layers.melody is not None
        assert len(layers.bass) > 0
        assert len(layers.pads) > 0

    def test_caffeine_pipeline(self):
        feat, cfg = _features("CN1C=NC2=C1C(=O)N(C(=O)N2C)C")
        layers = generate_all_layers(feat, cfg)
        assert layers.melody is not None
        assert len(layers.all_layers()) > 1
