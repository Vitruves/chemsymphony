"""§8: Functional group recognition via SMARTS patterns."""

from __future__ import annotations

from rdkit import Chem

from chemsymphony.features import MolecularFeatures

# (name, SMARTS pattern)
_FG_PATTERNS: list[tuple[str, str]] = [
    # Carboxylic acids & derivatives
    ("carboxylic_acid", "[CX3](=O)[OX2H1]"),
    ("ester", "[CX3](=O)[OX2][C]"),
    ("amide", "[CX3](=O)[NX3]"),
    ("acyl_halide", "[CX3](=O)[F,Cl,Br,I]"),
    ("anhydride", "[CX3](=O)[OX2][CX3](=O)"),
    ("carbamate", "[OX2][CX3](=O)[NX3]"),
    ("urea", "[NX3][CX3](=O)[NX3]"),
    ("lactone", "[C]1[C](=O)[O][$([C])]1"),
    ("lactam", "[C]1[C](=O)[N][$([C])]1"),
    # Carbonyls
    ("aldehyde", "[CX3H1](=O)[#6,H]"),
    ("ketone", "[#6][CX3](=O)[#6]"),
    # Alcohols & ethers
    ("hydroxyl", "[OX2H]"),
    ("phenol", "[OX2H][cX3]:[c]"),
    ("enol", "[OX2H][#6X3]=[#6]"),
    ("ether", "[OD2]([#6])[#6]"),
    ("epoxide", "[OX2r3]1[#6r3][#6r3]1"),
    ("peroxide", "[OX2][OX2]"),
    ("acetal", "[CX4]([OX2])([OX2])[H,#6]"),
    # Nitrogen-containing
    ("primary_amine", "[NX3H2][#6]"),
    ("secondary_amine", "[NX3H1]([#6])[#6]"),
    ("tertiary_amine", "[NX3]([#6])([#6])[#6]"),
    ("nitrile", "[CX2]#[NX1]"),
    ("isocyanate", "[NX2]=[CX2]=[OX1]"),
    ("isothiocyanate", "[NX2]=[CX2]=[SX1]"),
    ("nitro", "[$([NX3](=O)=O),$([NX3+](=O)[O-])]"),
    ("imine", "[CX3]([#6])=[NX2][#6,H]"),
    ("oxime", "[CX3]=[NX2][OX2H]"),
    ("hydrazine", "[NX3][NX3]"),
    ("hydrazone", "[CX3]=[NX2][NX3]"),
    ("azide", "[$([NX1]~[NX2]~[NX2]),$([NX2]=[NX2+]=[NX1-])]"),
    ("azo", "[#6][NX2]=[NX2][#6]"),
    ("guanidine", "[NX3][CX3](=[NX2])[NX3]"),
    ("enamine", "[NX3][CX3]=[CX3]"),
    # Sulfur-containing
    ("thiol", "[SX2H]"),
    ("thioether", "[SX2]([#6])[#6]"),
    ("disulfide", "[SX2][SX2]"),
    ("sulfoxide", "[SX3](=O)([#6])[#6]"),
    ("sulfone", "[SX4](=O)(=O)([#6])[#6]"),
    ("sulfonic_acid", "[SX4](=O)(=O)[OX2H]"),
    ("sulfonamide", "[SX4](=O)(=O)[NX3]"),
    ("sulfonate_ester", "[SX4](=O)(=O)[OX2][#6]"),
    ("thioamide", "[CX3](=S)[NX3]"),
    ("thioketone", "[#6][CX3](=S)[#6]"),
    # Phosphorus-containing
    ("phosphate", "[PX4](=O)([O])([O])[O]"),
    ("phosphonate", "[PX4](=O)([OX2])([OX2])[#6]"),
    ("phosphine", "[PX3]([#6])([#6])[#6]"),
    ("phosphoramide", "[PX4](=O)([NX3])([NX3])[NX3,OX2]"),
    # Boron-containing
    ("boronic_acid", "[BX3]([OX2H])[OX2H]"),
    ("boronate_ester", "[BX3]([OX2][#6])[OX2][#6]"),
    # Halides
    ("halide", "[#6][F,Cl,Br,I]"),
    # Unsaturated
    ("alkene", "[CX3]=[CX3]"),
    ("alkyne", "[CX2]#[CX2]"),
    # Heterocyclic
    ("pyridine", "[nX2]1[cX3][cX3][cX3][cX3][cX3]1"),
    ("pyrrole", "[nX3H1]1[cX3][cX3][cX3][cX3]1"),
    ("furan", "[oX2]1[cX3][cX3][cX3][cX3]1"),
    ("thiophene", "[sX2]1[cX3][cX3][cX3][cX3]1"),
    ("imidazole", "[nX3]1[cX3][nX2][cX3][cX3]1"),
    ("oxazole", "[oX2]1[cX3][nX2][cX3][cX3]1"),
    ("thiazole", "[sX2]1[cX3][nX2][cX3][cX3]1"),
    ("indole", "[nX3H1]1[cX3][cX3]2[cX3][cX3][cX3][cX3][cX3]2[cX3]1"),
]

# Pre-compile SMARTS
_COMPILED: list[tuple[str, Chem.Mol]] = []


def _get_compiled() -> list[tuple[str, Chem.Mol]]:
    if not _COMPILED:
        for name, smarts in _FG_PATTERNS:
            pattern = Chem.MolFromSmarts(smarts)
            if pattern is not None:
                _COMPILED.append((name, pattern))
    return _COMPILED


def extract_functional_groups(mol: Chem.Mol, feat: MolecularFeatures) -> None:
    """Populate functional group features on *feat*."""
    groups = []
    for name, pattern in _get_compiled():
        matches = mol.GetSubstructMatches(pattern)
        for match in matches:
            groups.append({
                "name": name,
                "atoms": list(match),
            })
    feat.functional_groups = groups
