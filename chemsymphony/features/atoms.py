"""§2: Atomic composition — element counts, fractions, unique elements."""

from __future__ import annotations

from collections import Counter

from rdkit import Chem

from chemsymphony.features import MolecularFeatures


def extract_atoms(mol: Chem.Mol, feat: MolecularFeatures) -> None:
    """Populate atomic composition features on *feat*."""
    symbols = [a.GetSymbol() for a in mol.GetAtoms()]
    counts = Counter(symbols)
    total = len(symbols) or 1

    feat.element_counts = dict(counts)
    feat.element_fractions = {el: c / total for el, c in counts.items()}
    feat.unique_elements = sorted(counts.keys())

    # Include implicit hydrogens in element_counts as well
    h_count = sum(a.GetTotalNumHs() for a in mol.GetAtoms())
    if h_count > 0:
        feat.element_counts["H"] = feat.element_counts.get("H", 0) + h_count
        total_with_h = total + h_count
        # Recompute fractions including H
        feat.element_fractions = {
            el: c / total_with_h for el, c in feat.element_counts.items()
        }
        if "H" not in feat.unique_elements:
            feat.unique_elements = sorted(feat.unique_elements + ["H"])
