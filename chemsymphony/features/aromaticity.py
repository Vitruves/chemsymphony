"""§4: Aromatic system analysis — aromatic ring count, conjugation, heteroaromatic."""

from __future__ import annotations

from rdkit import Chem

from chemsymphony.features import MolecularFeatures


def extract_aromaticity(mol: Chem.Mol, feat: MolecularFeatures) -> None:
    """Populate aromaticity-related features on *feat*."""
    aromatic_atoms = {a.GetIdx() for a in mol.GetAtoms() if a.GetIsAromatic()}
    feat.aromatic_atom_count = len(aromatic_atoms)

    total_heavy = mol.GetNumHeavyAtoms()
    feat.aromatic_fraction = len(aromatic_atoms) / total_heavy if total_heavy > 0 else 0.0

    # Count aromatic rings from the already-extracted ring data
    aromatic_ring_count = 0
    heteroaromatic_rings = []
    for ring in feat.rings:
        if ring["is_aromatic"]:
            aromatic_ring_count += 1
            heteroatoms_in_ring = ring["heteroatoms"]
            if heteroatoms_in_ring:
                heteroaromatic_rings.append({
                    "ring_index": ring["index"],
                    "size": ring["size"],
                    "heteroatoms": heteroatoms_in_ring,
                })

    feat.aromatic_ring_count = aromatic_ring_count
    feat.heteroaromatic_rings = heteroaromatic_rings

    # Conjugation length — longest path of aromatic atoms
    if not aromatic_atoms:
        feat.conjugation_length = 0
        return

    # BFS-based longest path within the aromatic subgraph
    adj: dict[int, list[int]] = {idx: [] for idx in aromatic_atoms}
    for bond in mol.GetBonds():
        a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if a in aromatic_atoms and b in aromatic_atoms:
            adj[a].append(b)
            adj[b].append(a)

    best = 0
    for start in aromatic_atoms:
        visited = {start}
        stack = [(start, 1)]
        while stack:
            node, depth = stack.pop()
            best = max(best, depth)
            for nb in adj[node]:
                if nb not in visited:
                    visited.add(nb)
                    stack.append((nb, depth + 1))

    feat.conjugation_length = best
