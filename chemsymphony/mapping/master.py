"""§1 master mapping: global molecular properties → audio parameters."""

from __future__ import annotations

from chemsymphony.config import Config
from chemsymphony.features import MolecularFeatures

# Scale definitions as semitone intervals from root
SCALES: dict[str, list[int]] = {
    "pentatonic_major": [0, 2, 4, 7, 9],
    "major": [0, 2, 4, 5, 7, 9, 11],
    "mixolydian": [0, 2, 4, 5, 7, 9, 10],
    "dorian": [0, 2, 3, 5, 7, 9, 10],
    "minor": [0, 2, 3, 5, 7, 8, 10],
    "phrygian": [0, 1, 3, 5, 7, 8, 10],
    "chromatic": [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11],
}

NOTE_NAMES = ["C", "C#", "D", "D#", "E", "F", "F#", "G", "G#", "A", "A#", "B"]


def _compute_bpm(mw: float, cfg: Config) -> int:
    """Map molecular weight to BPM."""
    if mw < 100:
        bpm = 70 + (mw / 100) * 20  # 70–90
    elif mw < 500:
        bpm = 90 + ((mw - 100) / 400) * 40  # 90–130
    else:
        bpm = 130 + min((mw - 500) / 500, 1.0) * 50  # 130–180
    return int(max(cfg.min_bpm, min(cfg.max_bpm, bpm)))


def _compute_duration(heavy_atom_count: int, cfg: Config) -> float:
    """Map heavy atom count to duration in seconds."""
    dur = 2.0 + 0.3 * heavy_atom_count
    return max(cfg.min_duration, min(cfg.max_duration, dur))


def _choose_scale(heteroatom_ratio: float) -> str:
    """Map heteroatom-to-carbon ratio to a musical scale."""
    if heteroatom_ratio <= 0.05:
        return "pentatonic_major"
    elif heteroatom_ratio <= 0.3:
        return "major"
    elif heteroatom_ratio <= 0.45:
        return "mixolydian"
    elif heteroatom_ratio <= 0.6:
        return "dorian"
    elif heteroatom_ratio <= 0.8:
        return "minor"
    elif heteroatom_ratio <= 1.0:
        return "phrygian"
    else:
        return "chromatic"


def _harmonic_density(unique_element_count: int) -> int:
    """Map unique element count to harmonic density (max simultaneous notes)."""
    if unique_element_count <= 2:
        return 1
    elif unique_element_count <= 4:
        return 3
    else:
        return min(unique_element_count, 6)


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
    bpm = cfg.bpm if cfg.bpm is not None else _compute_bpm(feat.molecular_weight, cfg)

    # Duration
    duration = cfg.duration if cfg.duration is not None else _compute_duration(
        feat.heavy_atom_count, cfg
    )

    # Root note and scale
    if cfg.key is not None:
        root_idx, forced_mode = _parse_key(cfg.key)
        scale_name = forced_mode or _choose_scale(feat.heteroatom_ratio)
    else:
        root_idx = feat.formula_hash  # 0–11
        scale_name = _choose_scale(feat.heteroatom_ratio)

    root_note_midi = 60 + root_idx  # Middle C octave
    scale_intervals = SCALES[scale_name]

    harmonic_dens = _harmonic_density(feat.unique_element_count)
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
