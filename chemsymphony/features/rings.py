"""§3: Ring detection and characterization — SSSR, fused/spiro, substituents."""

from __future__ import annotations

from rdkit import Chem

from chemsymphony.features import MolecularFeatures


def extract_rings(mol: Chem.Mol, feat: MolecularFeatures) -> None:
    """Populate ring-related features on *feat*."""
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()

    feat.ring_count = len(atom_rings)
    if feat.ring_count == 0:
        return

    rings = []
    for i, ring_atoms in enumerate(atom_rings):
        ring_atoms = list(ring_atoms)
        size = len(ring_atoms)

        # Atom composition of the ring
        symbols = [mol.GetAtomWithIdx(a).GetSymbol() for a in ring_atoms]
        heteroatoms = [s for s in symbols if s != "C"]
        is_aromatic = all(mol.GetAtomWithIdx(a).GetIsAromatic() for a in ring_atoms)

        # Count substituents: atoms bonded to ring atoms that are NOT in the ring
        ring_set = set(ring_atoms)
        substituent_count = 0
        for a_idx in ring_atoms:
            atom = mol.GetAtomWithIdx(a_idx)
            for neighbor in atom.GetNeighbors():
                if neighbor.GetIdx() not in ring_set:
                    substituent_count += 1

        rings.append({
            "index": i,
            "size": size,
            "atoms": ring_atoms,
            "symbols": symbols,
            "heteroatoms": heteroatoms,
            "is_aromatic": is_aromatic,
            "substituent_count": substituent_count,
        })

    feat.rings = rings

    # Fused ring detection (shared atoms between rings)
    fused_pairs = []
    spiro_atoms_set: set[int] = set()
    for i in range(len(atom_rings)):
        for j in range(i + 1, len(atom_rings)):
            shared = set(atom_rings[i]) & set(atom_rings[j])
            if len(shared) >= 2:
                fused_pairs.append((i, j))
            elif len(shared) == 1:
                spiro_atoms_set.update(shared)

    feat.fused_ring_pairs = fused_pairs
    feat.spiro_atoms = sorted(spiro_atoms_set)

    # Ring systems — group rings by connectivity
    # Build adjacency for rings
    ring_adj: dict[int, set[int]] = {i: set() for i in range(len(atom_rings))}
    for i in range(len(atom_rings)):
        for j in range(i + 1, len(atom_rings)):
            if set(atom_rings[i]) & set(atom_rings[j]):
                ring_adj[i].add(j)
                ring_adj[j].add(i)

    visited: set[int] = set()
    systems: list[list[int]] = []
    for start in range(len(atom_rings)):
        if start in visited:
            continue
        system = []
        stack = [start]
        while stack:
            node = stack.pop()
            if node in visited:
                continue
            visited.add(node)
            system.append(node)
            stack.extend(ring_adj[node] - visited)
        systems.append(sorted(system))

    feat.ring_systems = systems
