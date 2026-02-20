"""§7: Bond type analysis — counts, fractions, order sequence."""

from __future__ import annotations

from rdkit import Chem

from chemsymphony.features import MolecularFeatures


def extract_bonds(mol: Chem.Mol, feat: MolecularFeatures) -> None:
    """Populate bond-related features on *feat*."""
    single = double = triple = aromatic = 0

    for bond in mol.GetBonds():
        bt = bond.GetBondType()
        if bt == Chem.BondType.SINGLE:
            single += 1
        elif bt == Chem.BondType.DOUBLE:
            double += 1
        elif bt == Chem.BondType.TRIPLE:
            triple += 1
        elif bt == Chem.BondType.AROMATIC:
            aromatic += 1

    total = mol.GetNumBonds() or 1
    feat.single_bond_count = single
    feat.double_bond_count = double
    feat.triple_bond_count = triple
    feat.aromatic_bond_count = aromatic
    feat.double_bond_fraction = double / total
    feat.triple_bond_fraction = triple / total

    # Bond order sequence along the main chain (already computed in chains.py)
    # feat.bond_order_sequence comes from chain_bond_orders
    feat.bond_order_sequence = list(feat.chain_bond_orders)
