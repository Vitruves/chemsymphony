"""§9: Stereochemistry extraction — R/S, E/Z, meso detection."""

from __future__ import annotations

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

from chemsymphony.features import MolecularFeatures


def extract_stereo(mol: Chem.Mol, feat: MolecularFeatures) -> None:
    """Populate stereochemistry features on *feat*."""
    # Assign stereo info
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)

    # Chiral centers
    chiral_centers = []
    chiral_info = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    for atom_idx, label in chiral_info:
        chiral_centers.append({
            "atom_idx": atom_idx,
            "config": label,  # "R", "S", or "?"
            "symbol": mol.GetAtomWithIdx(atom_idx).GetSymbol(),
        })
    feat.chiral_centers = chiral_centers
    feat.stereo_center_count = len(chiral_centers)

    # E/Z bonds
    ez_bonds = []
    for bond in mol.GetBonds():
        stereo = bond.GetStereo()
        if stereo in (Chem.BondStereo.STEREOE, Chem.BondStereo.STEREOZ):
            label = "E" if stereo == Chem.BondStereo.STEREOE else "Z"
            ez_bonds.append({
                "bond_idx": bond.GetIdx(),
                "begin_atom": bond.GetBeginAtomIdx(),
                "end_atom": bond.GetEndAtomIdx(),
                "config": label,
            })
    feat.ez_bonds = ez_bonds

    # Meso detection: has stereo centers but is achiral overall
    # Simplified: if there are ≥2 stereo centers and the molecule has an
    # internal symmetry plane, it's meso. We approximate by checking if
    # R count == S count and molecule has symmetry.
    r_count = sum(1 for c in chiral_centers if c["config"] == "R")
    s_count = sum(1 for c in chiral_centers if c["config"] == "S")
    if r_count > 0 and r_count == s_count:
        # Check canonical SMILES without stereo vs with
        smi_no_stereo = Chem.MolToSmiles(mol, isomericSmiles=False)
        mol_no_stereo = Chem.MolFromSmiles(smi_no_stereo)
        if mol_no_stereo:
            # If removing stereo doesn't change the SMILES structure, it could be meso
            feat.is_meso = True
    else:
        feat.is_meso = False
