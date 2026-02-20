"""§1: Global molecular properties — heavy atom count, MW, heteroatom ratio, formula hash."""

from __future__ import annotations

import hashlib

from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors

from chemsymphony.features import MolecularFeatures


def extract_global_props(mol: Chem.Mol, feat: MolecularFeatures) -> None:
    """Populate global molecular properties on *feat*."""
    feat.heavy_atom_count = mol.GetNumHeavyAtoms()
    feat.molecular_weight = Descriptors.ExactMolWt(mol)

    formula = rdMolDescriptors.CalcMolFormula(mol)
    feat.molecular_formula = formula

    # Heteroatom-to-carbon ratio
    carbon_count = sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() == 6)
    heteroatom_count = feat.heavy_atom_count - carbon_count
    if carbon_count > 0:
        feat.heteroatom_ratio = heteroatom_count / carbon_count
    else:
        # No carbons — all heteroatoms
        feat.heteroatom_ratio = float(feat.heavy_atom_count) if feat.heavy_atom_count > 0 else 0.0

    # Unique element types (heavy atoms only)
    elements = {a.GetSymbol() for a in mol.GetAtoms()}
    feat.unique_element_count = len(elements)

    # Total bond count
    feat.total_bond_count = mol.GetNumBonds()

    # Deterministic formula hash → root note selection (0–11)
    h = hashlib.sha256(formula.encode("utf-8")).hexdigest()
    feat.formula_hash = int(h, 16) % 12
