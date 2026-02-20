"""Feature extraction pipelines for molecular analysis."""

from __future__ import annotations

from dataclasses import dataclass, field, fields
from typing import Any

from rdkit import Chem


@dataclass
class MolecularFeatures:
    """Aggregated features extracted from a molecule.

    Each field corresponds to the output of one feature-extraction pipeline.
    """

    # Global properties (§1)
    heavy_atom_count: int = 0
    molecular_weight: float = 0.0
    heteroatom_ratio: float = 0.0
    unique_element_count: int = 0
    total_bond_count: int = 0
    molecular_formula: str = ""
    formula_hash: int = 0

    # Atomic composition (§2)
    element_counts: dict[str, int] = field(default_factory=dict)
    element_fractions: dict[str, float] = field(default_factory=dict)
    unique_elements: list[str] = field(default_factory=list)

    # Rings (§3)
    ring_count: int = 0
    rings: list[dict[str, Any]] = field(default_factory=list)
    fused_ring_pairs: list[tuple[int, int]] = field(default_factory=list)
    spiro_atoms: list[int] = field(default_factory=list)
    ring_systems: list[list[int]] = field(default_factory=list)

    # Aromaticity (§4)
    aromatic_ring_count: int = 0
    aromatic_atom_count: int = 0
    aromatic_fraction: float = 0.0
    conjugation_length: int = 0
    heteroaromatic_rings: list[dict[str, Any]] = field(default_factory=list)

    # Chains (§5)
    longest_chain: int = 0
    longest_chain_atoms: list[int] = field(default_factory=list)
    branch_points_on_chain: list[int] = field(default_factory=list)
    chain_bond_orders: list[int] = field(default_factory=list)
    chain_heteroatom_positions: list[int] = field(default_factory=list)

    # Branches (§6)
    branch_count: int = 0
    branches: list[dict[str, Any]] = field(default_factory=list)

    # Bonds (§7)
    single_bond_count: int = 0
    double_bond_count: int = 0
    triple_bond_count: int = 0
    aromatic_bond_count: int = 0
    double_bond_fraction: float = 0.0
    triple_bond_fraction: float = 0.0
    bond_order_sequence: list[int] = field(default_factory=list)

    # Functional groups (§8)
    functional_groups: list[dict[str, Any]] = field(default_factory=list)

    # Stereochemistry (§9)
    chiral_centers: list[dict[str, Any]] = field(default_factory=list)
    ez_bonds: list[dict[str, Any]] = field(default_factory=list)
    is_meso: bool = False
    stereo_center_count: int = 0

    # Topology (§10)
    graph_diameter: int = 0
    avg_degree: float = 0.0
    max_degree: int = 0
    max_degree_atom: int = 0
    connected_components: int = 1
    symmetry_score: float = 0.0
    bridges: list[tuple[int, int]] = field(default_factory=list)
    wiener_index: int = 0

    # Distribution (§11)
    element_positions: dict[str, list[int]] = field(default_factory=dict)
    element_clustering: dict[str, int] = field(default_factory=dict)
    element_periodicity: dict[str, float] = field(default_factory=dict)
    element_ratios: dict[str, float] = field(default_factory=dict)

    # Electronic (§12)
    formal_charges: list[dict[str, Any]] = field(default_factory=dict)
    net_charge: int = 0
    is_zwitterion: bool = False
    electronegativity_gradient: list[float] = field(default_factory=list)
    radical_electrons: list[int] = field(default_factory=list)

    # Physicochemical (§13)
    logp: float = 0.0
    tpsa: float = 0.0
    rotatable_bond_count: int = 0
    hbd_count: int = 0
    hba_count: int = 0
    fsp3: float = 0.0
    bertz_ct: float = 0.0
    num_valence_electrons: int = 0
    num_radical_electrons_total: int = 0
    hall_kier_alpha: float = 0.0

    # Derived audio parameters (set by master mapping)
    audio_parameters: dict[str, Any] = field(default_factory=dict)

    def pretty(self) -> str:
        """Return a human-readable summary of all extracted features."""
        lines: list[str] = []
        a = lines.append

        def section(title: str) -> None:
            a("")
            a(f"--- {title} ---")

        def kv(key: str, val: object, indent: int = 2) -> None:
            a(f"{' ' * indent}{key}: {val}")

        # Header
        a(f"  {self.molecular_formula}  (MW {self.molecular_weight:.2f}, {self.heavy_atom_count} heavy atoms, {self.total_bond_count} bonds)")

        # Global
        section("Global Properties")
        kv("heavy_atom_count", self.heavy_atom_count)
        kv("molecular_weight", f"{self.molecular_weight:.4f}")
        kv("heteroatom_ratio", f"{self.heteroatom_ratio:.3f}")
        kv("unique_element_count", self.unique_element_count)
        kv("total_bond_count", self.total_bond_count)
        kv("molecular_formula", self.molecular_formula)
        kv("formula_hash", self.formula_hash)

        # Atoms
        section("Atomic Composition")
        kv("element_counts", self.element_counts)
        kv("element_fractions", {k: round(v, 4) for k, v in self.element_fractions.items()})
        kv("unique_elements", self.unique_elements)

        # Rings
        section("Rings")
        kv("ring_count", self.ring_count)
        for ring in self.rings:
            arom = "aromatic" if ring["is_aromatic"] else "saturated"
            hetero = f", hetero=[{', '.join(ring['heteroatoms'])}]" if ring["heteroatoms"] else ""
            a(f"    ring {ring['index']+1}: {ring['size']}-membered {arom}, subs={ring['substituent_count']}{hetero}")
        if self.fused_ring_pairs:
            kv("fused_ring_pairs", self.fused_ring_pairs)
        if self.spiro_atoms:
            kv("spiro_atoms", self.spiro_atoms)
        if self.ring_systems:
            kv("ring_systems", self.ring_systems)

        # Aromaticity
        section("Aromaticity")
        kv("aromatic_ring_count", self.aromatic_ring_count)
        kv("aromatic_atom_count", self.aromatic_atom_count)
        kv("aromatic_fraction", f"{self.aromatic_fraction:.3f}")
        kv("conjugation_length", self.conjugation_length)
        if self.heteroaromatic_rings:
            kv("heteroaromatic_rings", self.heteroaromatic_rings)

        # Chain
        section("Chain & Branches")
        kv("longest_chain", self.longest_chain)
        kv("longest_chain_atoms", self.longest_chain_atoms)
        if self.branch_points_on_chain:
            kv("branch_points_on_chain", self.branch_points_on_chain)
        kv("chain_bond_orders", self.chain_bond_orders)
        if self.chain_heteroatom_positions:
            kv("chain_heteroatom_positions", self.chain_heteroatom_positions)
        kv("branch_count", self.branch_count)
        for b in self.branches:
            sym_tag = " (symmetric)" if b.get("is_symmetric") else ""
            a(f"    branch @{b['chain_position']}: len={b['length']}, depth={b['depth']}, atoms={b['symbols']}{sym_tag}")

        # Bonds
        section("Bonds")
        kv("single", self.single_bond_count)
        kv("double", f"{self.double_bond_count}  ({self.double_bond_fraction:.3f})")
        kv("triple", f"{self.triple_bond_count}  ({self.triple_bond_fraction:.3f})")
        kv("aromatic", self.aromatic_bond_count)
        kv("bond_order_sequence", self.bond_order_sequence)

        # Functional groups
        section("Functional Groups")
        if self.functional_groups:
            seen: dict[str, int] = {}
            for fg in self.functional_groups:
                seen[fg["name"]] = seen.get(fg["name"], 0) + 1
            for name, count in seen.items():
                suffix = f" x{count}" if count > 1 else ""
                a(f"    {name}{suffix}")
        else:
            a("    (none)")

        # Stereo
        section("Stereochemistry")
        kv("stereo_center_count", self.stereo_center_count)
        for c in self.chiral_centers:
            a(f"    {c['symbol']}@{c['atom_idx']}: {c['config']}")
        for ez in self.ez_bonds:
            a(f"    bond {ez['begin_atom']}-{ez['end_atom']}: {ez['config']}")
        kv("is_meso", self.is_meso)

        # Topology
        section("Topology")
        kv("graph_diameter", self.graph_diameter)
        kv("avg_degree", f"{self.avg_degree:.3f}")
        kv("max_degree", f"{self.max_degree} (atom {self.max_degree_atom})")
        kv("connected_components", self.connected_components)
        kv("symmetry_score", f"{self.symmetry_score:.3f}")
        kv("wiener_index", self.wiener_index)
        kv("bridges", self.bridges)

        # Distribution
        section("Distribution")
        kv("element_positions", self.element_positions)
        kv("element_clustering", self.element_clustering)
        kv("element_periodicity", {k: round(v, 3) for k, v in self.element_periodicity.items()})
        kv("element_ratios", {k: round(v, 3) for k, v in self.element_ratios.items()})

        # Electronic
        section("Electronic")
        if self.formal_charges:
            kv("formal_charges", self.formal_charges)
        kv("net_charge", self.net_charge)
        kv("is_zwitterion", self.is_zwitterion)
        if self.electronegativity_gradient:
            kv("electronegativity_gradient", [round(v, 2) for v in self.electronegativity_gradient])
        if self.radical_electrons:
            kv("radical_electrons", self.radical_electrons)

        # Physicochemical
        section("Physicochemical Properties")
        kv("logp", f"{self.logp:.3f}")
        kv("tpsa", f"{self.tpsa:.2f}")
        kv("rotatable_bond_count", self.rotatable_bond_count)
        kv("hbd_count", self.hbd_count)
        kv("hba_count", self.hba_count)
        kv("fsp3", f"{self.fsp3:.3f}")
        kv("bertz_ct", f"{self.bertz_ct:.2f}")
        kv("num_valence_electrons", self.num_valence_electrons)
        kv("num_radical_electrons_total", self.num_radical_electrons_total)
        kv("hall_kier_alpha", f"{self.hall_kier_alpha:.3f}")

        # Audio parameters
        ap = self.audio_parameters
        if ap:
            section("Audio Parameters")
            kv("bpm", ap.get("bpm"))
            kv("duration", f"{ap.get('duration', 0):.1f}s")
            kv("root_note", f"{ap.get('root_note')} (MIDI {ap.get('root_note_midi')})")
            kv("scale", ap.get("scale"))
            kv("scale_intervals", ap.get("scale_intervals"))
            kv("harmonic_density", ap.get("harmonic_density"))
            kv("note_density", f"{ap.get('note_density', 0):.2f}")
            kv("filter_warmth", f"{ap.get('filter_warmth', 0):.3f}")
            kv("reverb_wetness", f"{ap.get('reverb_wetness', 0):.3f}")
            kv("swing", f"{ap.get('swing', 0):.3f}")
            kv("timbre_organic", f"{ap.get('timbre_organic', 0):.3f}")
            kv("arrangement_density", f"{ap.get('arrangement_density', 1.0):.3f}")
            kv("harmonic_tension", f"{ap.get('harmonic_tension', 0):.3f}")
            if ap.get("seed") is not None:
                kv("seed", ap["seed"])

        return "\n".join(lines)


def extract_all_features(mol: Chem.Mol) -> MolecularFeatures:
    """Run all 12 feature extraction pipelines and return aggregated results."""
    from chemsymphony.features.global_props import extract_global_props
    from chemsymphony.features.atoms import extract_atoms
    from chemsymphony.features.rings import extract_rings
    from chemsymphony.features.aromaticity import extract_aromaticity
    from chemsymphony.features.chains import extract_chains
    from chemsymphony.features.branches import extract_branches
    from chemsymphony.features.bonds import extract_bonds
    from chemsymphony.features.functional_groups import extract_functional_groups
    from chemsymphony.features.stereo import extract_stereo
    from chemsymphony.features.topology import extract_topology
    from chemsymphony.features.distribution import extract_distribution
    from chemsymphony.features.electronic import extract_electronic
    from chemsymphony.features.physicochemical import extract_physicochemical

    feat = MolecularFeatures()

    extract_global_props(mol, feat)
    extract_atoms(mol, feat)
    extract_rings(mol, feat)
    extract_aromaticity(mol, feat)
    extract_chains(mol, feat)
    extract_branches(mol, feat)
    extract_bonds(mol, feat)
    extract_functional_groups(mol, feat)
    extract_stereo(mol, feat)
    extract_topology(mol, feat)
    extract_distribution(mol, feat)
    extract_electronic(mol, feat)
    extract_physicochemical(mol, feat)

    return feat
