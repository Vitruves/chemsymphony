"""§15 feature extraction: Simple graph counting features.

Light-weight features computed by iterating atoms and bonds without
expensive RDKit descriptor calculations.
"""

from __future__ import annotations

import itertools
import statistics
from collections import Counter

from rdkit import Chem
from rdkit.Chem import HybridizationType

from chemsymphony.features import MolecularFeatures

# Pauling electronegativity values
_EN: dict[str, float] = {
    "H": 2.20, "C": 2.55, "N": 3.04, "O": 3.44, "F": 3.98,
    "P": 2.19, "S": 2.58, "Cl": 3.16, "Se": 2.55, "Br": 2.96,
    "I": 2.66, "Si": 1.90, "B": 2.04, "Na": 0.93, "K": 0.82,
    "Mg": 1.31, "Ca": 1.00, "Fe": 1.83, "Zn": 1.65, "Cu": 1.90,
}


def extract_graph_counting(mol: Chem.Mol, feat: MolecularFeatures) -> None:
    """Extract simple counting features from the molecular graph."""

    atoms = list(mol.GetAtoms())

    # 1. Atom hybridization histogram
    feat.sp_count = sum(
        1 for a in atoms if a.GetHybridization() == HybridizationType.SP
    )
    feat.sp2_count = sum(
        1 for a in atoms if a.GetHybridization() == HybridizationType.SP2
    )
    feat.sp3_count = sum(
        1 for a in atoms if a.GetHybridization() == HybridizationType.SP3
    )

    # 2. Heteroatom adjacency count (bonds where both atoms are non-carbon)
    feat.heteroatom_adjacency_count = sum(
        1 for b in mol.GetBonds()
        if b.GetBeginAtom().GetAtomicNum() != 6
        and b.GetEndAtom().GetAtomicNum() != 6
    )

    # 3. Terminal atom count (degree-1 heavy atoms)
    feat.terminal_atom_count = sum(1 for a in atoms if a.GetDegree() == 1)

    # 4. Quaternary center count (atoms bonded to 4+ heavy atoms)
    feat.quaternary_center_count = sum(
        1 for a in atoms
        if sum(1 for n in a.GetNeighbors() if n.GetAtomicNum() > 1) >= 4
    )

    # 5. Neighbor element pair frequencies
    pair_counts: Counter[str] = Counter()
    for b in mol.GetBonds():
        s1 = b.GetBeginAtom().GetSymbol()
        s2 = b.GetEndAtom().GetSymbol()
        pair = "-".join(sorted([s1, s2]))
        pair_counts[pair] += 1
    feat.neighbor_pair_counts = dict(pair_counts)

    # 6. Chain-to-ring atom ratio
    ring_info = mol.GetRingInfo()
    ring_atom_set: set[int] = set()
    for ring in ring_info.AtomRings():
        ring_atom_set.update(ring)
    ring_atom_count = len(ring_atom_set)
    chain_atom_count = feat.heavy_atom_count - ring_atom_count
    if ring_atom_count > 0:
        feat.chain_to_ring_ratio = chain_atom_count / ring_atom_count
    else:
        feat.chain_to_ring_ratio = float(chain_atom_count) if chain_atom_count > 0 else 0.0

    # 7. Maximum consecutive same-element run along longest chain
    if feat.longest_chain_atoms:
        symbols = [
            mol.GetAtomWithIdx(idx).GetSymbol()
            for idx in feat.longest_chain_atoms
        ]
        feat.max_same_element_run = max(
            (len(list(g)) for _, g in itertools.groupby(symbols)), default=0
        )
    else:
        feat.max_same_element_run = 0

    # 8. Functional group diversity
    if feat.functional_groups:
        fg_names = {fg["name"] for fg in feat.functional_groups}
        feat.fg_diversity = len(fg_names)
    else:
        feat.fg_diversity = 0

    # 9. Ring size variance
    if feat.rings and len(feat.rings) > 1:
        sizes = [r["size"] for r in feat.rings]
        feat.ring_size_variance = statistics.variance(sizes)
    else:
        feat.ring_size_variance = 0.0

    # 10. Atom electronegativity variance
    en_values = [
        _EN.get(a.GetSymbol(), 2.2) for a in atoms
    ]
    if len(en_values) > 1:
        feat.en_variance = statistics.variance(en_values)
    else:
        feat.en_variance = 0.0

    # 11. Macrocycle detection (ring with >12 atoms)
    feat.has_macrocycle = any(r["size"] > 12 for r in feat.rings) if feat.rings else False
