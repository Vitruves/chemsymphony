"""§12: Charges & electronegativity — formal charges, net charge, zwitterion, radicals."""

from __future__ import annotations

from rdkit import Chem

from chemsymphony.features import MolecularFeatures

# Pauling electronegativity values for common elements
_ELECTRONEGATIVITY: dict[int, float] = {
    1: 2.20,   # H
    6: 2.55,   # C
    7: 3.04,   # N
    8: 3.44,   # O
    9: 3.98,   # F
    15: 2.19,  # P
    16: 2.58,  # S
    17: 3.16,  # Cl
    35: 2.96,  # Br
    53: 2.66,  # I
    11: 0.93,  # Na
    19: 0.82,  # K
    26: 1.83,  # Fe
}


def extract_electronic(mol: Chem.Mol, feat: MolecularFeatures) -> None:
    """Populate electronic features on *feat*."""
    charges = []
    has_positive = False
    has_negative = False
    net = 0
    radical_positions = []

    for atom in mol.GetAtoms():
        fc = atom.GetFormalCharge()
        if fc != 0:
            charges.append({
                "atom_idx": atom.GetIdx(),
                "symbol": atom.GetSymbol(),
                "charge": fc,
            })
            net += fc
            if fc > 0:
                has_positive = True
            else:
                has_negative = True

        rad = atom.GetNumRadicalElectrons()
        if rad > 0:
            radical_positions.append(atom.GetIdx())

    feat.formal_charges = charges
    feat.net_charge = net
    feat.is_zwitterion = has_positive and has_negative
    feat.radical_electrons = radical_positions

    # Electronegativity gradient along the longest chain
    gradient = []
    for idx in feat.longest_chain_atoms:
        atom = mol.GetAtomWithIdx(idx)
        en = _ELECTRONEGATIVITY.get(atom.GetAtomicNum(), 2.5)
        gradient.append(en)
    feat.electronegativity_gradient = gradient
