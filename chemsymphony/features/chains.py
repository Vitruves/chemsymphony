"""§5: Longest chain extraction — BFS longest path, branch points, bond orders."""

from __future__ import annotations

from collections import deque

from rdkit import Chem

from chemsymphony.features import MolecularFeatures


def _bfs_farthest(adj: dict[int, list[int]], start: int) -> tuple[int, dict[int, int]]:
    """BFS from *start*, return (farthest_node, distance_map)."""
    dist: dict[int, int] = {start: 0}
    queue = deque([start])
    farthest = start
    while queue:
        node = queue.popleft()
        for nb in adj[node]:
            if nb not in dist:
                dist[nb] = dist[node] + 1
                queue.append(nb)
                if dist[nb] > dist[farthest]:
                    farthest = nb
    return farthest, dist


def _reconstruct_path(adj: dict[int, list[int]], start: int, end: int) -> list[int]:
    """BFS shortest path from *start* to *end*."""
    if start == end:
        return [start]
    parent: dict[int, int] = {start: -1}
    queue = deque([start])
    while queue:
        node = queue.popleft()
        for nb in adj[node]:
            if nb not in parent:
                parent[nb] = node
                if nb == end:
                    path = []
                    cur = end
                    while cur != -1:
                        path.append(cur)
                        cur = parent[cur]
                    return path[::-1]
                queue.append(nb)
    return [start]


def extract_chains(mol: Chem.Mol, feat: MolecularFeatures) -> None:
    """Populate longest-chain features on *feat*."""
    n = mol.GetNumAtoms()
    if n == 0:
        return

    # Build adjacency list (all heavy atoms)
    adj: dict[int, list[int]] = {i: [] for i in range(n)}
    bond_order_map: dict[tuple[int, int], int] = {}
    for bond in mol.GetBonds():
        a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        adj[a].append(b)
        adj[b].append(a)
        bt = bond.GetBondType()
        order = {
            Chem.BondType.SINGLE: 1,
            Chem.BondType.DOUBLE: 2,
            Chem.BondType.TRIPLE: 3,
            Chem.BondType.AROMATIC: 1,
        }.get(bt, 1)
        bond_order_map[(a, b)] = order
        bond_order_map[(b, a)] = order

    if n == 1:
        feat.longest_chain = 1
        feat.longest_chain_atoms = [0]
        return

    # Double-BFS to find diameter path (longest shortest path)
    far1, _ = _bfs_farthest(adj, 0)
    far2, _ = _bfs_farthest(adj, far1)
    path = _reconstruct_path(adj, far1, far2)

    feat.longest_chain = len(path)
    feat.longest_chain_atoms = path

    # Branch points on the chain: chain atoms with degree > 2 in the full graph
    chain_set = set(path)
    branch_points = []
    for idx in path:
        non_chain_neighbors = sum(1 for nb in adj[idx] if nb not in chain_set)
        if non_chain_neighbors > 0 and len(adj[idx]) > 2:
            branch_points.append(idx)
    feat.branch_points_on_chain = branch_points

    # Bond orders along the chain
    bond_orders = []
    for i in range(len(path) - 1):
        order = bond_order_map.get((path[i], path[i + 1]), 1)
        bond_orders.append(order)
    feat.chain_bond_orders = bond_orders

    # Heteroatom positions in chain
    hetero_positions = []
    for i, idx in enumerate(path):
        atom = mol.GetAtomWithIdx(idx)
        if atom.GetAtomicNum() != 6:
            hetero_positions.append(i)
    feat.chain_heteroatom_positions = hetero_positions
