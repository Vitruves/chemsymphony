"""§13 feature extraction: Physicochemical properties via RDKit descriptors."""

from __future__ import annotations

from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors, GraphDescriptors


def extract_physicochemical(mol: Chem.Mol, feat: object) -> None:
    """Extract physicochemical properties and store on *feat*."""
    feat.logp = Descriptors.MolLogP(mol)
    feat.tpsa = Descriptors.TPSA(mol)
    feat.rotatable_bond_count = rdMolDescriptors.CalcNumRotatableBonds(mol)
    feat.hbd_count = rdMolDescriptors.CalcNumHBD(mol)
    feat.hba_count = rdMolDescriptors.CalcNumHBA(mol)
    feat.fsp3 = rdMolDescriptors.CalcFractionCSP3(mol)
    feat.bertz_ct = GraphDescriptors.BertzCT(mol)
    feat.num_valence_electrons = Descriptors.NumValenceElectrons(mol)

    # Total radical electrons across all atoms
    feat.num_radical_electrons_total = sum(
        a.GetNumRadicalElectrons() for a in mol.GetAtoms()
    )

    feat.hall_kier_alpha = GraphDescriptors.HallKierAlpha(mol)
