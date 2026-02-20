"""§10: Graph topology analysis — diameter, degrees, bridges, Wiener index, symmetry."""

from __future__ import annotations

from collections import deque

from rdkit import Chem

from chemsymphony.features import MolecularFeatures


def _bfs_distances(adj: dict[int, list[int]], start: int) -> dict[int, int]:
    dist = {start: 0}
    queue = deque([start])
    while queue:
        node = queue.popleft()
        for nb in adj[node]:
            if nb not in dist:
                dist[nb] = dist[node] + 1
                queue.append(nb)
    return dist


def _find_bridges(adj: dict[int, list[int]], n: int) -> list[tuple[int, int]]:
    """Find bridges using Tarjan's algorithm."""
    disc: dict[int, int] = {}
    low: dict[int, int] = {}
    bridges: list[tuple[int, int]] = []
    timer = [0]

    def dfs(u: int, parent: int) -> None:
        disc[u] = low[u] = timer[0]
        timer[0] += 1
        for v in adj[u]:
            if v not in disc:
                dfs(v, u)
                low[u] = min(low[u], low[v])
                if low[v] > disc[u]:
                    bridges.append((u, v))
            elif v != parent:
                low[u] = min(low[u], disc[v])

    for node in adj:
        if node not in disc:
            dfs(node, -1)

    return bridges


def extract_topology(mol: Chem.Mol, feat: MolecularFeatures) -> None:
    """Populate graph topology features on *feat*."""
    n = mol.GetNumAtoms()
    if n == 0:
        return

    adj: dict[int, list[int]] = {i: [] for i in range(n)}
    for bond in mol.GetBonds():
        a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        adj[a].append(b)
        adj[b].append(a)

    # Degrees
    degrees = [len(adj[i]) for i in range(n)]
    feat.avg_degree = sum(degrees) / n if n > 0 else 0.0
    feat.max_degree = max(degrees) if degrees else 0
    feat.max_degree_atom = degrees.index(feat.max_degree) if degrees else 0

    # Connected components
    visited: set[int] = set()
    components = 0
    for i in range(n):
        if i not in visited:
            components += 1
            queue = deque([i])
            visited.add(i)
            while queue:
                node = queue.popleft()
                for nb in adj[node]:
                    if nb not in visited:
                        visited.add(nb)
                        queue.append(nb)
    feat.connected_components = components

    # BFS from all nodes to compute diameter and Wiener index
    diameter = 0
    wiener = 0
    for i in range(n):
        dists = _bfs_distances(adj, i)
        if dists:
            max_dist = max(dists.values())
            diameter = max(diameter, max_dist)
            wiener += sum(dists.values())
    feat.graph_diameter = diameter
    feat.wiener_index = wiener // 2  # Each pair counted twice

    # Bridges (Tarjan's)
    feat.bridges = _find_bridges(adj, n)

    # Symmetry score: approximate via canonical SMILES atom map
    # Use RDKit's canonical ranking as a proxy
    ranks = list(Chem.CanonicalRankAtoms(mol))
    unique_ranks = len(set(ranks))
    feat.symmetry_score = 1.0 - (unique_ranks / n) if n > 1 else 0.0
