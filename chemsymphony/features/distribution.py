"""§11: Element distribution — positions in traversal, clustering, periodicity, ratios."""

from __future__ import annotations

from collections import Counter, defaultdict

from rdkit import Chem

from chemsymphony.features import MolecularFeatures


def extract_distribution(mol: Chem.Mol, feat: MolecularFeatures) -> None:
    """Populate element distribution features on *feat*."""
    n = mol.GetNumAtoms()
    if n == 0:
        return

    # Element positions in canonical traversal order (atom index order)
    positions: dict[str, list[int]] = defaultdict(list)
    for idx in range(n):
        sym = mol.GetAtomWithIdx(idx).GetSymbol()
        positions[sym].append(idx)
    feat.element_positions = dict(positions)

    # Element clustering — longest run of adjacent same-element atoms
    adj: dict[int, list[int]] = {i: [] for i in range(n)}
    for bond in mol.GetBonds():
        a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        adj[a].append(b)
        adj[b].append(a)

    clustering: dict[str, int] = {}
    for sym, idxs in positions.items():
        if len(idxs) <= 1:
            clustering[sym] = 1
            continue
        idx_set = set(idxs)
        # BFS within same-element subgraph to find largest component
        visited: set[int] = set()
        max_cluster = 0
        for start in idxs:
            if start in visited:
                continue
            cluster_size = 0
            queue = [start]
            visited.add(start)
            while queue:
                node = queue.pop()
                cluster_size += 1
                for nb in adj[node]:
                    if nb in idx_set and nb not in visited:
                        visited.add(nb)
                        queue.append(nb)
            max_cluster = max(max_cluster, cluster_size)
        clustering[sym] = max_cluster
    feat.element_clustering = clustering

    # Element periodicity — regularity of spacing
    periodicity: dict[str, float] = {}
    for sym, idxs in positions.items():
        if len(idxs) < 2:
            periodicity[sym] = 0.0
            continue
        gaps = [idxs[i + 1] - idxs[i] for i in range(len(idxs) - 1)]
        mean_gap = sum(gaps) / len(gaps)
        if mean_gap == 0:
            periodicity[sym] = 1.0
        else:
            variance = sum((g - mean_gap) ** 2 for g in gaps) / len(gaps)
            # Regularity = 1 / (1 + coefficient_of_variation)
            std = variance ** 0.5
            periodicity[sym] = 1.0 / (1.0 + std / mean_gap) if mean_gap > 0 else 0.0
    feat.element_periodicity = periodicity

    # Element ratios (relative to the most common heavy element)
    counts = Counter(mol.GetAtomWithIdx(i).GetSymbol() for i in range(n))
    if counts:
        max_count = max(counts.values())
        feat.element_ratios = {sym: c / max_count for sym, c in counts.items()}
    else:
        feat.element_ratios = {}
