"""§14 feature extraction: SMILES string-level features.

These features operate on the canonical SMILES text representation,
capturing representational complexity without molecular graph analysis.
"""

from __future__ import annotations

import math
import re
import statistics
from collections import Counter

from chemsymphony.features import MolecularFeatures


def extract_smiles_features(smiles: str, feat: MolecularFeatures) -> None:
    """Extract features from the canonical SMILES string itself."""
    feat.canonical_smiles = smiles

    if not smiles:
        return

    # 1. Shannon entropy of character distribution
    counts = Counter(smiles)
    n = len(smiles)
    feat.smiles_entropy = -sum(
        (c / n) * math.log2(c / n) for c in counts.values()
    ) if n > 0 else 0.0

    # 2. Max parenthesis nesting depth
    depth = max_d = 0
    for ch in smiles:
        if ch == '(':
            depth += 1
            max_d = max(max_d, depth)
        elif ch == ')':
            depth -= 1
    feat.smiles_nesting_depth = max_d

    # 3. SMILES length / heavy atom count ratio
    hac = max(feat.heavy_atom_count, 1)
    feat.smiles_length_ratio = len(smiles) / hac

    # 4. Bracket atom count — atoms requiring explicit [...]
    feat.bracket_atom_count = len(re.findall(r'\[.*?\]', smiles))

    # 5. Maximum ring closure number
    singles = [int(d) for d in re.findall(r'(?<![%\d])(\d)(?!\d)', smiles)]
    multis = [int(m) for m in re.findall(r'%(\d+)', smiles)]
    all_nums = singles + multis
    feat.max_ring_closure = max(all_nums) if all_nums else 0

    # 6. Aromatic character ratio (lowercase c/n/o/s vs uppercase C/N/O/S)
    lower_aromatic = sum(1 for ch in smiles if ch in 'cnos')
    upper_aliphatic = sum(1 for ch in smiles if ch in 'CNOS')
    total_atoms_chars = lower_aromatic + upper_aliphatic
    feat.aromatic_char_ratio = (
        lower_aromatic / total_atoms_chars if total_atoms_chars > 0 else 0.0
    )

    # 7. Repeating substrings of length >= 3
    best_motif = ""
    best_count = 0
    total_repeats = 0
    seen: set[str] = set()
    for length in range(3, min(len(smiles) // 2 + 1, 12)):
        for i in range(len(smiles) - length + 1):
            sub = smiles[i:i + length]
            if sub in seen:
                continue
            seen.add(sub)
            count = smiles.count(sub)
            if count > 1:
                total_repeats += 1
                if len(sub) * count > len(best_motif) * best_count:
                    best_motif = sub
                    best_count = count
    feat.repeating_motif_count = total_repeats
    feat.longest_repeating_motif = best_motif

    # 8. Special character density
    special = sum(1 for ch in smiles if ch in '=#@/\\+-[]()%.')
    feat.special_char_density = special / n if n > 0 else 0.0

    # 9. Fragment length variance (dot-separated fragments)
    fragments = smiles.split('.')
    if len(fragments) > 1:
        lengths = [len(f) for f in fragments]
        feat.fragment_length_variance = statistics.variance(lengths)
    else:
        feat.fragment_length_variance = 0.0
