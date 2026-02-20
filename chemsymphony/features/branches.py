"""§6: Branch detection — substituent trees off the longest chain."""

from __future__ import annotations

from collections import deque

from rdkit import Chem

from chemsymphony.features import MolecularFeatures


def extract_branches(mol: Chem.Mol, feat: MolecularFeatures) -> None:
    """Populate branch-related features on *feat*."""
    chain_atoms = set(feat.longest_chain_atoms)
    if not chain_atoms:
        return

    n = mol.GetNumAtoms()
    adj: dict[int, list[int]] = {i: [] for i in range(n)}
    for bond in mol.GetBonds():
        a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        adj[a].append(b)
        adj[b].append(a)

    branches = []

    for chain_pos, chain_idx in enumerate(feat.longest_chain_atoms):
        for nb in adj[chain_idx]:
            if nb in chain_atoms:
                continue
            # This is the root of a branch
            branch_atoms = []
            visited = set(chain_atoms)
            queue = deque([nb])
            visited.add(nb)
            max_depth = 0
            depth_map = {nb: 1}

            while queue:
                node = queue.popleft()
                branch_atoms.append(node)
                for n2 in adj[node]:
                    if n2 not in visited:
                        visited.add(n2)
                        depth_map[n2] = depth_map[node] + 1
                        max_depth = max(max_depth, depth_map[n2])
                        queue.append(n2)

            symbols = [mol.GetAtomWithIdx(a).GetSymbol() for a in branch_atoms]
            branches.append({
                "root_atom": nb,
                "chain_position": chain_pos,
                "chain_atom": chain_idx,
                "length": len(branch_atoms),
                "depth": max_depth if max_depth else 1,
                "atoms": branch_atoms,
                "symbols": symbols,
            })

    # Detect symmetric (identical) branches by their symbol tuples
    symbol_tuples = [tuple(sorted(b["symbols"])) for b in branches]
    for i, branch in enumerate(branches):
        branch["is_symmetric"] = symbol_tuples.count(symbol_tuples[i]) > 1

    feat.branch_count = len(branches)
    feat.branches = branches
