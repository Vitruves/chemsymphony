"""Tests for CLI."""

import tempfile
from pathlib import Path

import pytest
from chemsymphony.cli import main


class TestCLI:
    def test_help(self, capsys):
        with pytest.raises(SystemExit) as exc_info:
            main(["--help"])
        assert exc_info.value.code == 0
        captured = capsys.readouterr()
        assert "SMILES" in captured.out or "smiles" in captured.out

    def test_version(self, capsys):
        with pytest.raises(SystemExit) as exc_info:
            main(["--version"])
        assert exc_info.value.code == 0
        captured = capsys.readouterr()
        assert "0.3.0" in captured.out

    def test_dry_run(self, capsys):
        main(["-s", "CCO", "--dry-run"])
        captured = capsys.readouterr()
        assert "heavy_atom_count" in captured.out
        assert "Audio Parameters" in captured.out

    def test_invalid_smiles(self):
        with pytest.raises(SystemExit) as exc_info:
            main(["-s", "INVALID_SMILES"])
        assert exc_info.value.code != 0

    def test_generate_wav(self):
        with tempfile.NamedTemporaryFile(suffix=".wav", delete=False) as f:
            path = f.name
        main(["-s", "CCO", "-f", "wav", "-o", path])
        assert Path(path).stat().st_size > 44
        Path(path).unlink()

    def test_generate_midi(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            outdir = str(Path(tmpdir) / "test_midi")
            main(["-s", "c1ccccc1", "-f", "midi", "-o", outdir])
            assert (Path(outdir) / "master.mid").exists()

    def test_verbose(self, capsys):
        with tempfile.NamedTemporaryFile(suffix=".wav", delete=False) as f:
            path = f.name
        main(["-s", "CCO", "-f", "wav", "-o", path, "-v"])
        captured = capsys.readouterr()
        assert "heavy_atom_count" in captured.out
        Path(path).unlink()

    def test_quiet(self, capsys):
        with tempfile.NamedTemporaryFile(suffix=".wav", delete=False) as f:
            path = f.name
        main(["-s", "CCO", "-f", "wav", "-o", path, "-q"])
        captured = capsys.readouterr()
        assert captured.out.strip() == ""
        Path(path).unlink()

    def test_bpm_override(self, capsys):
        main(["-s", "CCO", "--dry-run", "--bpm", "120"])

    def test_seed(self, capsys):
        main(["-s", "CCO", "--dry-run", "--seed", "42"])

    def test_missing_smiles(self):
        with pytest.raises(SystemExit) as exc_info:
            main([])
        assert exc_info.value.code != 0
