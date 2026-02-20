"""Tests for export functionality."""

import json
import tempfile
from pathlib import Path

import pytest
from chemsymphony.core import ChemSymphony, MidiResult
from chemsymphony.config import load_config


class TestMidiExport:
    def test_generate_midi(self):
        cs = ChemSymphony()
        result = cs.generate("c1ccccc1", fmt="midi")
        assert isinstance(result, MidiResult)
        assert len(result.master) > 0
        assert len(result.tracks) > 0
        assert "smiles_canonical" in result.manifest

    def test_save_midi_directory(self):
        cs = ChemSymphony()
        result = cs.generate("c1ccccc1", fmt="midi")
        with tempfile.TemporaryDirectory() as tmpdir:
            outdir = Path(tmpdir) / "benzene"
            result.save(outdir)
            assert (outdir / "master.mid").exists()
            assert (outdir / "manifest.json").exists()
            assert (outdir / "README.txt").exists()
            assert (outdir / "tracks").is_dir()

            # Verify manifest is valid JSON
            manifest = json.loads((outdir / "manifest.json").read_text())
            assert manifest["smiles_canonical"] == "c1ccccc1"


class TestWavExport:
    def test_generate_wav(self):
        cs = ChemSymphony()
        result = cs.generate("CCO", fmt="wav")
        assert isinstance(result, bytes)
        assert len(result) > 44  # WAV header is 44 bytes
        assert result[:4] == b"RIFF"
        assert result[8:12] == b"WAVE"

    def test_generate_to_file(self):
        cs = ChemSymphony()
        with tempfile.NamedTemporaryFile(suffix=".wav", delete=False) as f:
            path = f.name
        cs.generate_to_file("CCO", output=path, fmt="wav")
        data = Path(path).read_bytes()
        assert data[:4] == b"RIFF"
        Path(path).unlink()


class TestDryRun:
    def test_extract_features(self):
        cs = ChemSymphony()
        feat = cs.extract_features("CC(=O)OC1=CC=CC=C1C(=O)O")
        assert feat.ring_count >= 1
        assert feat.heavy_atom_count > 0
        assert "bpm" in feat.audio_parameters
        assert feat.molecular_formula != ""
