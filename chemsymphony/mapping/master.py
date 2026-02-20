"""§1 master mapping: global molecular properties → audio parameters."""

from __future__ import annotations

import math

from chemsymphony.config import Config
from chemsymphony.features import MolecularFeatures

# Scale definitions as semitone intervals from root
SCALES: dict[str, list[int]] = {
    "pentatonic_major": [0, 2, 4, 7, 9],
    "pentatonic_minor": [0, 3, 5, 7, 10],
    "major": [0, 2, 4, 5, 7, 9, 11],
    "lydian": [0, 2, 4, 6, 7, 9, 11],
    "mixolydian": [0, 2, 4, 5, 7, 9, 10],
    "dorian": [0, 2, 3, 5, 7, 9, 10],
    "melodic_minor": [0, 2, 3, 5, 7, 9, 11],
    "minor": [0, 2, 3, 5, 7, 8, 10],
    "harmonic_minor": [0, 2, 3, 5, 7, 8, 11],
    "blues": [0, 3, 5, 6, 7, 10],
    "phrygian": [0, 1, 3, 5, 7, 8, 10],
    "whole_tone": [0, 2, 4, 6, 8, 10],
    "hungarian_minor": [0, 2, 3, 6, 7, 8, 11],
    "chromatic": [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11],
    "locrian": [0, 1, 3, 5, 6, 8, 10],
    "lydian_dominant": [0, 2, 4, 6, 7, 9, 10],
    "altered": [0, 1, 3, 4, 6, 8, 10],
    "harmonic_major": [0, 2, 4, 5, 7, 8, 11],
    "double_harmonic": [0, 1, 4, 5, 7, 8, 11],
    "bebop_dominant": [0, 2, 4, 5, 7, 9, 10, 11],
}

NOTE_NAMES = ["C", "C#", "D", "D#", "E", "F", "F#", "G", "G#", "A", "A#", "B"]


def _sigmoid(x: float, center: float, steepness: float = 1.0) -> float:
    """Sigmoid mapping to (0, 1), centered at *center*."""
    return 1.0 / (1.0 + math.exp(-steepness * (x - center)))


def _compute_bpm(feat: MolecularFeatures, cfg: Config) -> int:
    """Map molecular properties to BPM using 5-axis weighted scoring.

    Axes (weights):
      MW (25%)             — sigmoid centred at 300 for wider spread
      aromatic_fraction (20%) — conjugation → driving rhythm
      bertz_ct (20%)       — complexity → faster tempo (log-scaled)
      graph_diameter (20%, inv) — elongation → slower
      symmetry_score (15%) — high symmetry → mid-tempo
    Secondary modulations: rotatable_bonds, fsp3, ring_density (±10 BPM).
    """
    # --- primary axes (each normalised to 0–1) ---
    # MW: sigmoid centred at 300, steepness 0.008
    mw_score = _sigmoid(feat.molecular_weight, 300.0, 0.008)

    # Aromatic fraction (already 0–1)
    arom_score = feat.aromatic_fraction

    # Bertz complexity (log-scaled, ~100–2000 range)
    bertz = max(feat.bertz_ct, 1.0)
    bertz_score = min(1.0, math.log(bertz) / math.log(2000))

    # Graph diameter (inverted — elongated = slower)
    diameter_score = 1.0 - min(1.0, feat.graph_diameter / 20.0)

    # Symmetry (maps to mid-tempo: 0.5 at extremes, 1.0 at symmetry=0.5)
    sym = feat.symmetry_score
    sym_score = 1.0 - 2.0 * abs(sym - 0.5)  # Peak at 0.5

    # Weighted sum → 0–1
    composite = (
        0.25 * mw_score
        + 0.20 * arom_score
        + 0.20 * bertz_score
        + 0.20 * diameter_score
        + 0.15 * sym_score
    )

    # Map composite to BPM range (60–180)
    bpm = 60 + composite * 120  # 60–180

    # --- secondary modulations (±10 BPM) ---
    bpm += min(feat.rotatable_bond_count * 0.8, 5.0)
    bpm -= 6 * feat.fsp3
    if feat.heavy_atom_count > 0:
        ring_density = feat.ring_count / feat.heavy_atom_count
        bpm += 8 * ring_density

    return int(max(cfg.min_bpm, min(cfg.max_bpm, bpm)))


def _compute_duration(heavy_atom_count: int, cfg: Config) -> float:
    """Map heavy atom count to duration in seconds."""
    dur = 2.0 + 0.3 * heavy_atom_count
    return max(cfg.min_duration, min(cfg.max_duration, dur))


def _choose_scale(heteroatom_ratio: float,
                   aromatic_fraction: float = 0.0,
                   fsp3: float = 0.0,
                   ring_count: int = 0,
                   fg_diversity: int = 0) -> str:
    """Map heteroatom ratio + aromatic fraction + fsp3 to a musical scale.

    Uses a 2D decision space: heteroatom_ratio controls darkness/tension,
    aromatic_fraction and fsp3 select between related scales.
    ring_count and fg_diversity (unique FG types) further subdivide buckets.
    """
    if heteroatom_ratio <= 0.05:
        return "pentatonic_major"
    elif heteroatom_ratio <= 0.15:
        if fsp3 > 0.5:
            return "pentatonic_minor"
        return "major"
    elif heteroatom_ratio <= 0.25:
        if aromatic_fraction > 0.4 and ring_count >= 3 and fg_diversity >= 3:
            return "bebop_dominant"
        if aromatic_fraction > 0.4:
            return "mixolydian"
        elif fsp3 > 0.3:
            return "melodic_minor"
        return "lydian"
    elif heteroatom_ratio <= 0.4:
        if aromatic_fraction > 0.5 and fg_diversity >= 3:
            return "altered"
        if aromatic_fraction > 0.5:
            return "dorian"
        elif ring_count >= 3 and aromatic_fraction > 0.2:
            return "lydian_dominant"
        elif fsp3 > 0.4:
            return "blues"
        return "mixolydian"
    elif heteroatom_ratio <= 0.6:
        if aromatic_fraction > 0.3 and fg_diversity >= 4:
            return "altered"
        if aromatic_fraction > 0.3:
            return "harmonic_minor"
        if ring_count >= 2:
            return "harmonic_major"
        return "dorian"
    elif heteroatom_ratio <= 0.8:
        if fsp3 > 0.5:
            return "whole_tone"
        if fg_diversity >= 3:
            return "locrian"
        return "minor"
    elif heteroatom_ratio <= 1.0:
        if aromatic_fraction > 0.3:
            return "hungarian_minor"
        if fg_diversity >= 3:
            return "double_harmonic"
        return "phrygian"
    else:
        return "chromatic"


def _harmonic_density(unique_element_count: int,
                      fused_ring_pair_count: int = 0,
                      branch_count: int = 0) -> int:
    """Map unique elements + fused rings + branches to harmonic density."""
    base = unique_element_count
    # Fused rings add harmonic complexity
    base += fused_ring_pair_count
    # Branches add polyphonic texture
    base += branch_count // 2

    if base <= 2:
        return 1
    elif base <= 3:
        return 2
    elif base <= 5:
        return 3
    elif base <= 7:
        return 4
    elif base <= 9:
        return 5
    else:
        return 6


def _note_density(total_bonds: int, duration: float) -> float:
    """Map total bond count to events per beat."""
    # Base: 1 event/beat, scale up with bonds
    return max(0.5, min(4.0, total_bonds / max(duration, 1.0)))


def _parse_key(key_str: str) -> tuple[int, str | None]:
    """Parse a key string like 'Cm', 'G', 'F#m' into (root_midi_offset, mode).

    Returns (root_note_index_0_11, mode_or_None).
    """
    key_str = key_str.strip()
    mode = None
    if key_str.endswith("m"):
        mode = "minor"
        key_str = key_str[:-1]

    note_map = {
        "C": 0, "C#": 1, "Db": 1, "D": 2, "D#": 3, "Eb": 3,
        "E": 4, "F": 5, "F#": 6, "Gb": 6, "G": 7, "G#": 8,
        "Ab": 8, "A": 9, "A#": 10, "Bb": 10, "B": 11,
    }
    root = note_map.get(key_str, 0)
    return root, mode


def _clamp(val: float, lo: float, hi: float) -> float:
    return max(lo, min(hi, val))


def apply_master_mapping(feat: MolecularFeatures, cfg: Config) -> None:
    """Set audio_parameters on *feat* based on global molecular properties."""
    # BPM
    bpm = cfg.bpm if cfg.bpm is not None else _compute_bpm(feat, cfg)

    # Duration
    duration = cfg.duration if cfg.duration is not None else _compute_duration(
        feat.heavy_atom_count, cfg
    )

    # Compute FG diversity (unique functional group types)
    fg_names = {fg["name"] for fg in feat.functional_groups} if feat.functional_groups else set()
    fg_diversity = len(fg_names)

    # Root note and scale
    if cfg.key is not None:
        root_idx, forced_mode = _parse_key(cfg.key)
        scale_name = forced_mode or _choose_scale(
            feat.heteroatom_ratio, feat.aromatic_fraction, feat.fsp3,
            ring_count=feat.ring_count, fg_diversity=fg_diversity,
        )
    else:
        root_idx = feat.formula_hash  # 0–11
        scale_name = _choose_scale(
            feat.heteroatom_ratio, feat.aromatic_fraction, feat.fsp3,
            ring_count=feat.ring_count, fg_diversity=fg_diversity,
        )

    root_note_midi = 60 + root_idx  # Middle C octave
    scale_intervals = SCALES[scale_name]

    harmonic_dens = _harmonic_density(
        feat.unique_element_count,
        len(feat.fused_ring_pairs),
        feat.branch_count,
    )
    note_dens = _note_density(feat.total_bond_count, duration)

    # §13 physicochemical → audio parameter derivations
    filter_warmth = _clamp(feat.logp / 5.0, 0.0, 1.0)
    reverb_wetness = _clamp(feat.tpsa / 140.0, 0.0, 1.0)
    swing = _clamp(feat.rotatable_bond_count * 0.03, 0.0, 0.3)
    timbre_organic = _clamp(feat.fsp3, 0.0, 1.0)
    arrangement_density = _clamp(0.5 + feat.bertz_ct / 2000.0, 0.5, 1.5)
    harmonic_tension = _clamp((feat.hbd_count + feat.hba_count) / 20.0, 0.0, 1.0)

    # §14 SMILES string → audio parameter derivations
    timbral_richness = _clamp(feat.smiles_entropy / 4.0, 0.0, 1.0)
    reverb_diffusion = _clamp(feat.smiles_nesting_depth / 5.0, 0.0, 1.0)
    ornament_density = _clamp(feat.special_char_density / 0.4, 0.0, 1.0)
    analog_warmth = _clamp(feat.bracket_atom_count / 5.0, 0.0, 1.0)

    # §15 graph counting → audio parameter derivations
    total_hyb = max(feat.sp_count + feat.sp2_count + feat.sp3_count, 1)
    waveform_brightness = _clamp(
        (feat.sp_count * 1.0 + feat.sp2_count * 0.5) / total_hyb, 0.0, 1.0
    )
    chromatic_dissonance = _clamp(feat.heteroatom_adjacency_count / 8.0, 0.0, 1.0)
    echo_density = _clamp(feat.terminal_atom_count / 10.0, 0.0, 1.0)
    filter_sweep = _clamp(feat.en_variance / 0.3, 0.0, 1.0)

    feat.audio_parameters = {
        "bpm": bpm,
        "duration": duration,
        "root_note": NOTE_NAMES[root_idx],
        "root_note_midi": root_note_midi,
        "scale": scale_name,
        "scale_intervals": scale_intervals,
        "harmonic_density": harmonic_dens,
        "note_density": note_dens,
        "seed": cfg.seed,
        "filter_warmth": filter_warmth,
        "reverb_wetness": reverb_wetness,
        "swing": swing,
        "timbre_organic": timbre_organic,
        "arrangement_density": arrangement_density,
        "harmonic_tension": harmonic_tension,
        "_symmetry_score": feat.symmetry_score,
        # §14 SMILES-derived
        "timbral_richness": timbral_richness,
        "reverb_diffusion": reverb_diffusion,
        "ornament_density": ornament_density,
        "analog_warmth": analog_warmth,
        # §15 graph-counting-derived
        "waveform_brightness": waveform_brightness,
        "chromatic_dissonance": chromatic_dissonance,
        "echo_density": echo_density,
        "filter_sweep": filter_sweep,
    }
