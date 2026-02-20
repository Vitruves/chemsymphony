"""Tests for SMILES canonicalization."""

import pytest
from chemsymphony.canonicalize import canonicalize, InvalidSmilesError


class TestCanonicalize:
    def test_simple_smiles(self):
        result = canonicalize("CCO")
        assert result.canonical_smiles == "CCO"
        assert result.mol is not None

    def test_equivalent_smiles_produce_same_canonical(self):
        r1 = canonicalize("OCC")
        r2 = canonicalize("C(O)C")
        r3 = canonicalize("CCO")
        assert r1.canonical_smiles == r2.canonical_smiles == r3.canonical_smiles

    def test_benzene(self):
        result = canonicalize("c1ccccc1")
        assert result.canonical_smiles == "c1ccccc1"

    def test_caffeine(self):
        result = canonicalize("CN1C=NC2=C1C(=O)N(C(=O)N2C)C")
        assert result.mol is not None
        assert len(result.canonical_smiles) > 0

    def test_invalid_smiles_raises(self):
        with pytest.raises(InvalidSmilesError):
            canonicalize("not_a_smiles")

    def test_empty_smiles_raises(self):
        with pytest.raises(InvalidSmilesError):
            canonicalize("")

    def test_whitespace_only_raises(self):
        with pytest.raises(InvalidSmilesError):
            canonicalize("   ")

    def test_single_atom(self):
        result = canonicalize("C")
        assert result.canonical_smiles == "C"

    def test_preserves_input(self):
        result = canonicalize("OCC")
        assert result.input_smiles == "OCC"
        assert result.canonical_smiles == "CCO"

    def test_water(self):
        result = canonicalize("O")
        assert result.canonical_smiles == "O"

    def test_salt(self):
        result = canonicalize("[Na+].[Cl-]")
        assert result.mol is not None
