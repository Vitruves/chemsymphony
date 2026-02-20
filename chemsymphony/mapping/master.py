"""§1 master mapping: global molecular properties → audio parameters."""

from __future__ import annotations

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
}

NOTE_NAMES = ["C", "C#", "D", "D#", "E", "F", "F#", "G", "G#", "A", "A#", "B"]


def _compute_bpm(feat: MolecularFeatures, cfg: Config) -> int:
    """Map molecular weight + structural features to BPM.

    Base BPM from MW, then modulated by rotatable bonds (+flexibility),
    fsp3 (-saturation slows tempo), and ring density (+compactness).
    """
    mw = feat.molecular_weight
    if mw < 100:
        bpm = 70 + (mw / 100) * 20  # 70–90
    elif mw < 500:
        bpm = 90 + ((mw - 100) / 400) * 40  # 90–130
    else:
        bpm = 130 + min((mw - 500) / 500, 1.0) * 50  # 130–180

    # Rotatable bonds add rhythmic motion
    bpm += feat.rotatable_bond_count * 1.5
    # High sp3 fraction → more saturated/organic → slower feel
    bpm -= 12 * feat.fsp3
    # Ring density → compact molecules feel faster
    if feat.heavy_atom_count > 0:
        ring_density = feat.ring_count / feat.heavy_atom_count
        bpm += 20 * ring_density

    return int(max(cfg.min_bpm, min(cfg.max_bpm, bpm)))


def _compute_duration(heavy_atom_count: int, cfg: Config) -> float:
    """Map heavy atom count to duration in seconds."""
    dur = 2.0 + 0.3 * heavy_atom_count
    return max(cfg.min_duration, min(cfg.max_duration, dur))


def _choose_scale(heteroatom_ratio: float,
                   aromatic_fraction: float = 0.0,
                   fsp3: float = 0.0) -> str:
    """Map heteroatom ratio + aromatic fraction + fsp3 to a musical scale.

    Uses a 2D decision space: heteroatom_ratio controls darkness/tension,
    aromatic_fraction and fsp3 select between related scales.
    """
    if heteroatom_ratio <= 0.05:
        return "pentatonic_major"
    elif heteroatom_ratio <= 0.15:
        # Low heteroatom: bright scales
        if fsp3 > 0.5:
            return "pentatonic_minor"  # Saturated + few heteroatoms
        return "major"
    elif heteroatom_ratio <= 0.25:
        # Moderate-low heteroatom
        if aromatic_fraction > 0.4:
            return "mixolydian"  # Aromatic dominant
        elif fsp3 > 0.3:
            return "melodic_minor"  # Mixed sp3 + heteroatoms
        return "lydian"
    elif heteroatom_ratio <= 0.4:
        if aromatic_fraction > 0.5:
            return "dorian"
        elif fsp3 > 0.4:
            return "blues"
        return "mixolydian"
    elif heteroatom_ratio <= 0.6:
        if aromatic_fraction > 0.3:
            return "harmonic_minor"
        return "dorian"
    elif heteroatom_ratio <= 0.8:
        if fsp3 > 0.5:
            return "whole_tone"
        return "minor"
    elif heteroatom_ratio <= 1.0:
        if aromatic_fraction > 0.3:
            return "hungarian_minor"
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

    # Root note and scale
    if cfg.key is not None:
        root_idx, forced_mode = _parse_key(cfg.key)
        scale_name = forced_mode or _choose_scale(
            feat.heteroatom_ratio, feat.aromatic_fraction, feat.fsp3
        )
    else:
        root_idx = feat.formula_hash  # 0–11
        scale_name = _choose_scale(
            feat.heteroatom_ratio, feat.aromatic_fraction, feat.fsp3
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
    }
