"""Tests for feature extraction pipelines."""

import pytest
from rdkit import Chem
from chemsymphony.features import extract_all_features


def _mol(smiles: str) -> Chem.Mol:
    return Chem.MolFromSmiles(smiles)


class TestGlobalProps:
    def test_ethanol(self):
        feat = extract_all_features(_mol("CCO"))
        assert feat.heavy_atom_count == 3
        assert feat.molecular_weight > 40
        assert feat.molecular_formula == "C2H6O"
        assert feat.unique_element_count == 2  # C, O (heavy atoms)
        assert 0 < feat.heteroatom_ratio < 1.0

    def test_methane(self):
        feat = extract_all_features(_mol("C"))
        assert feat.heavy_atom_count == 1
        assert feat.heteroatom_ratio == 0.0

    def test_water(self):
        feat = extract_all_features(_mol("O"))
        assert feat.heavy_atom_count == 1
        assert feat.heteroatom_ratio > 0

    def test_formula_hash_deterministic(self):
        f1 = extract_all_features(_mol("CCO"))
        f2 = extract_all_features(_mol("CCO"))
        assert f1.formula_hash == f2.formula_hash


class TestAtoms:
    def test_element_counts(self):
        feat = extract_all_features(_mol("CCO"))
        assert feat.element_counts["C"] == 2
        assert feat.element_counts["O"] == 1
        assert "H" in feat.element_counts

    def test_unique_elements(self):
        feat = extract_all_features(_mol("CCO"))
        assert "C" in feat.unique_elements
        assert "O" in feat.unique_elements


class TestRings:
    def test_benzene(self):
        feat = extract_all_features(_mol("c1ccccc1"))
        assert feat.ring_count == 1
        assert feat.rings[0]["size"] == 6
        assert feat.rings[0]["is_aromatic"] is True

    def test_naphthalene(self):
        feat = extract_all_features(_mol("c1ccc2ccccc2c1"))
        assert feat.ring_count == 2
        assert len(feat.fused_ring_pairs) >= 1

    def test_no_rings(self):
        feat = extract_all_features(_mol("CCCC"))
        assert feat.ring_count == 0

    def test_cyclohexane(self):
        feat = extract_all_features(_mol("C1CCCCC1"))
        assert feat.ring_count == 1
        assert feat.rings[0]["is_aromatic"] is False


class TestAromaticity:
    def test_benzene_aromatic(self):
        feat = extract_all_features(_mol("c1ccccc1"))
        assert feat.aromatic_ring_count == 1
        assert feat.aromatic_atom_count == 6
        assert feat.aromatic_fraction > 0.9

    def test_cyclohexane_not_aromatic(self):
        feat = extract_all_features(_mol("C1CCCCC1"))
        assert feat.aromatic_ring_count == 0

    def test_pyridine_heteroaromatic(self):
        feat = extract_all_features(_mol("c1ccncc1"))
        assert len(feat.heteroaromatic_rings) == 1


class TestChains:
    def test_butane(self):
        feat = extract_all_features(_mol("CCCC"))
        assert feat.longest_chain == 4

    def test_single_atom(self):
        feat = extract_all_features(_mol("C"))
        assert feat.longest_chain == 1

    def test_ethanol_chain(self):
        feat = extract_all_features(_mol("CCO"))
        assert feat.longest_chain >= 2


class TestBonds:
    def test_ethane(self):
        feat = extract_all_features(_mol("CC"))
        assert feat.single_bond_count >= 1
        assert feat.double_bond_count == 0

    def test_ethylene(self):
        feat = extract_all_features(_mol("C=C"))
        assert feat.double_bond_count >= 1

    def test_acetylene(self):
        feat = extract_all_features(_mol("C#C"))
        assert feat.triple_bond_count >= 1


class TestFunctionalGroups:
    def test_ethanol_hydroxyl(self):
        feat = extract_all_features(_mol("CCO"))
        names = [fg["name"] for fg in feat.functional_groups]
        assert "hydroxyl" in names

    def test_acetic_acid(self):
        feat = extract_all_features(_mol("CC(=O)O"))
        names = [fg["name"] for fg in feat.functional_groups]
        assert "carboxylic_acid" in names or "hydroxyl" in names


class TestStereo:
    def test_no_stereo(self):
        feat = extract_all_features(_mol("CC"))
        assert feat.stereo_center_count == 0

    def test_alanine(self):
        feat = extract_all_features(_mol("[C@@H](N)(C)C(=O)O"))
        assert feat.stereo_center_count >= 1


class TestTopology:
    def test_chain_diameter(self):
        feat = extract_all_features(_mol("CCCCCCCC"))
        assert feat.graph_diameter >= 7

    def test_single_atom_topology(self):
        feat = extract_all_features(_mol("C"))
        assert feat.graph_diameter == 0

    def test_salt_components(self):
        feat = extract_all_features(_mol("[Na+].[Cl-]"))
        assert feat.connected_components == 2


class TestPhysicochemical:
    """Tests for §13 physicochemical property extraction."""

    def test_ethanol_logp(self):
        feat = extract_all_features(_mol("CCO"))
        # Ethanol is hydrophilic, LogP should be negative/low
        assert feat.logp < 1.0

    def test_benzene_logp(self):
        feat = extract_all_features(_mol("c1ccccc1"))
        # Benzene is hydrophobic, LogP > 1
        assert feat.logp > 1.0

    def test_ethanol_tpsa(self):
        feat = extract_all_features(_mol("CCO"))
        # Ethanol has an OH group, so TPSA > 0
        assert feat.tpsa > 0

    def test_hexane_tpsa(self):
        feat = extract_all_features(_mol("CCCCCC"))
        # Pure hydrocarbon: TPSA = 0
        assert feat.tpsa == 0.0

    def test_rotatable_bonds(self):
        feat = extract_all_features(_mol("CCCCCC"))
        # Hexane has rotatable bonds
        assert feat.rotatable_bond_count >= 2

    def test_methane_no_rotatable(self):
        feat = extract_all_features(_mol("C"))
        assert feat.rotatable_bond_count == 0

    def test_fsp3_hexane(self):
        feat = extract_all_features(_mol("CCCCCC"))
        # All sp3 carbons
        assert feat.fsp3 == 1.0

    def test_fsp3_benzene(self):
        feat = extract_all_features(_mol("c1ccccc1"))
        # All sp2 carbons, no sp3
        assert feat.fsp3 == 0.0

    def test_bertz_ct_positive(self):
        feat = extract_all_features(_mol("CCO"))
        # Should be a positive complexity measure
        assert feat.bertz_ct > 0

    def test_hbd_hba_ethanol(self):
        feat = extract_all_features(_mol("CCO"))
        # Ethanol: 1 HBD (OH), 1 HBA (O)
        assert feat.hbd_count >= 1
        assert feat.hba_count >= 1

    def test_valence_electrons(self):
        feat = extract_all_features(_mol("CCO"))
        # C2H6O: C(4)*2 + H(1)*6 + O(6) = 20
        assert feat.num_valence_electrons == 20

    def test_radical_electrons_zero(self):
        feat = extract_all_features(_mol("CCO"))
        assert feat.num_radical_electrons_total == 0

    def test_hall_kier_alpha(self):
        feat = extract_all_features(_mol("CCO"))
        # Should return a numeric value
        assert isinstance(feat.hall_kier_alpha, float)

    def test_caffeine_properties(self):
        """Caffeine: moderate complexity molecule."""
        feat = extract_all_features(_mol("CN1C=NC2=C1C(=O)N(C(=O)N2C)C"))
        assert feat.logp > -2.0
        assert feat.tpsa > 30  # Multiple N and O atoms
        assert feat.bertz_ct > 100  # Complex molecule
        assert feat.hba_count >= 3  # Multiple N/O acceptors
