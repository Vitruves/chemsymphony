"""§8: Functional group recognition via SMARTS patterns."""

from __future__ import annotations

from rdkit import Chem

from chemsymphony.features import MolecularFeatures

# (name, SMARTS pattern)
_FG_PATTERNS: list[tuple[str, str]] = [
    ("carboxylic_acid", "[CX3](=O)[OX2H1]"),
    ("ester", "[CX3](=O)[OX2][C]"),
    ("amide", "[CX3](=O)[NX3]"),
    ("aldehyde", "[CX3H1](=O)[#6,H]"),
    ("ketone", "[#6][CX3](=O)[#6]"),
    ("hydroxyl", "[OX2H]"),
    ("ether", "[OD2]([#6])[#6]"),
    ("primary_amine", "[NX3H2][#6]"),
    ("secondary_amine", "[NX3H1]([#6])[#6]"),
    ("tertiary_amine", "[NX3]([#6])([#6])[#6]"),
    ("nitrile", "[CX2]#[NX1]"),
    ("nitro", "[$([NX3](=O)=O),$([NX3+](=O)[O-])]"),
    ("thiol", "[SX2H]"),
    ("sulfoxide", "[SX3](=O)([#6])[#6]"),
    ("phosphate", "[PX4](=O)([O])([O])[O]"),
    ("halide", "[#6][F,Cl,Br,I]"),
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
