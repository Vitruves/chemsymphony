"""Genre post-processing system.

Provides genre-specific audio styling via a scannable directory of genre modules.
Each module exports a ``GENRE = GenreProfile(...)`` object.  New genres are added
by dropping a Python file into the ``genres/`` package.

Two-phase approach:
  - **Pre-process** — clamp BPM range and adjust mix levels *before* generation
  - **Post-process** — apply EQ, compression, saturation, filtering, reverb, delay
    *after* ``mix_layers()`` returns the final stereo array
"""

from __future__ import annotations

import importlib
import pkgutil
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
from scipy import signal as scipy_signal

from chemsymphony.synthesis.effects import (
    apply_delay,
    apply_eq,
    apply_reverb,
    apply_soft_distortion,
    fade_in_out,
)


# ---------------------------------------------------------------------------
# GenreProfile dataclass
# ---------------------------------------------------------------------------

@dataclass
class GenreProfile:
    """Describes the sonic characteristics of a musical genre."""

    name: str

    # Pre-generation adjustments
    bpm_range: tuple[int, int] | None = None
    scale_preferences: list[str] | None = None
    mix_adjustments: dict[str, float] = field(default_factory=dict)
    swing_override: float | None = None

    # Post-processing EQ
    eq: dict[str, float] = field(default_factory=lambda: {
        "low_boost": 0.0, "mid_boost": 0.0, "brightness": 0.0,
    })

    # Post-processing dynamics
    compression: dict[str, float] = field(default_factory=lambda: {
        "threshold_db": -12.0, "ratio": 4.0,
    })

    # Post-processing effects
    reverb_depth: float | None = None
    saturation: float = 0.0
    highpass_hz: float | None = None
    lowpass_hz: float | None = None
    delay: dict[str, float] | None = None


# ---------------------------------------------------------------------------
# Genre discovery
# ---------------------------------------------------------------------------

def list_genres() -> list[str]:
    """Scan the ``genres/`` package for modules that export a ``GENRE`` attribute.

    Returns a sorted list of genre names (the module file names, not the
    ``GenreProfile.name`` display names).
    """
    package_path = str(Path(__file__).parent)
    names: list[str] = []
    for importer, modname, ispkg in pkgutil.iter_modules([package_path]):
        if ispkg or modname.startswith("_"):
            continue
        try:
            mod = importlib.import_module(f"chemsymphony.genres.{modname}")
            if hasattr(mod, "GENRE"):
                names.append(modname)
        except Exception:
            continue
    return sorted(names)


def load_genre(name: str) -> GenreProfile:
    """Import a genre module by name and return its ``GENRE`` object."""
    mod = importlib.import_module(f"chemsymphony.genres.{name}")
    profile = getattr(mod, "GENRE", None)
    if profile is None:
        raise ValueError(f"Genre module {name!r} does not export a GENRE object")
    return profile


# ---------------------------------------------------------------------------
# Pre-processing (applied before audio generation)
# ---------------------------------------------------------------------------

def apply_genre_preprocess(profile: GenreProfile, cfg: "Config") -> None:  # noqa: F821
    """Clamp BPM range, adjust mix levels, and override swing."""
    from chemsymphony.config import Config  # avoid circular import

    # Clamp BPM range
    if profile.bpm_range is not None:
        lo, hi = profile.bpm_range
        cfg.min_bpm = max(cfg.min_bpm, lo)
        cfg.max_bpm = min(cfg.max_bpm, hi)
        # If user set an explicit BPM, clamp it too
        if cfg.bpm is not None:
            cfg.bpm = max(lo, min(hi, cfg.bpm))

    # Multiply mix levels by genre adjustments
    for layer_name, multiplier in profile.mix_adjustments.items():
        if layer_name in cfg.mix:
            cfg.mix[layer_name] = cfg.mix[layer_name] * multiplier

    # Override swing
    if profile.swing_override is not None:
        cfg._genre_swing_override = profile.swing_override


# ---------------------------------------------------------------------------
# Post-processing (applied after mix_layers)
# ---------------------------------------------------------------------------

def _highpass(audio: np.ndarray, sr: int, cutoff_hz: float) -> np.ndarray:
    """Apply a Butterworth high-pass filter."""
    nyquist = sr / 2
    norm_cutoff = cutoff_hz / nyquist
    if norm_cutoff <= 0 or norm_cutoff >= 1:
        return audio
    sos = scipy_signal.butter(2, norm_cutoff, btype="high", output="sos")
    return scipy_signal.sosfilt(sos, audio)


def _lowpass(audio: np.ndarray, sr: int, cutoff_hz: float) -> np.ndarray:
    """Apply a Butterworth low-pass filter."""
    nyquist = sr / 2
    norm_cutoff = cutoff_hz / nyquist
    if norm_cutoff <= 0 or norm_cutoff >= 1:
        return audio
    sos = scipy_signal.butter(2, norm_cutoff, btype="low", output="sos")
    return scipy_signal.sosfilt(sos, audio)


def _compress_channel(audio: np.ndarray, sr: int,
                      threshold_db: float, ratio: float) -> np.ndarray:
    """Simple RMS compressor for a mono channel."""
    from chemsymphony.synthesis.mixer import _compress
    return _compress(audio, sr, threshold_db=threshold_db, ratio=ratio)


def _soft_limit_channel(audio: np.ndarray) -> np.ndarray:
    """Soft-limit a mono channel."""
    from chemsymphony.synthesis.mixer import _soft_limit
    return _soft_limit(audio)


def apply_genre_postprocess(
    audio: np.ndarray,
    profile: GenreProfile,
    sr: int,
    bpm: int,
) -> np.ndarray:
    """Apply genre-specific post-processing to a stereo array (N, 2).

    Chain: highpass -> lowpass -> EQ -> compression -> saturation ->
           reverb -> delay -> normalize -> soft-limit -> fade
    """
    left = audio[:, 0].copy()
    right = audio[:, 1].copy()

    # 1. High-pass filter
    if profile.highpass_hz is not None and profile.highpass_hz > 0:
        left = _highpass(left, sr, profile.highpass_hz)
        right = _highpass(right, sr, profile.highpass_hz)

    # 2. Low-pass filter
    if profile.lowpass_hz is not None and profile.lowpass_hz > 0:
        left = _lowpass(left, sr, profile.lowpass_hz)
        right = _lowpass(right, sr, profile.lowpass_hz)

    # 3. EQ
    eq = profile.eq
    if any(abs(v) >= 0.05 for v in eq.values()):
        left = apply_eq(left, sr,
                        brightness=eq.get("brightness", 0.0),
                        low_boost=eq.get("low_boost", 0.0),
                        mid_boost=eq.get("mid_boost", 0.0))
        right = apply_eq(right, sr,
                         brightness=eq.get("brightness", 0.0),
                         low_boost=eq.get("low_boost", 0.0),
                         mid_boost=eq.get("mid_boost", 0.0))

    # 4. Compression
    comp = profile.compression
    threshold = comp.get("threshold_db", -12.0)
    ratio = comp.get("ratio", 4.0)
    if ratio > 1.0:
        left = _compress_channel(left, sr, threshold, ratio)
        right = _compress_channel(right, sr, threshold, ratio)

    # 5. Saturation
    if profile.saturation > 0:
        left = apply_soft_distortion(left, amount=profile.saturation)
        right = apply_soft_distortion(right, amount=profile.saturation)

    # 6. Reverb override
    if profile.reverb_depth is not None and profile.reverb_depth > 0:
        left = apply_reverb(left, sr, depth=profile.reverb_depth)
        right = apply_reverb(right, sr, depth=profile.reverb_depth)

    # 7. Delay
    if profile.delay is not None:
        feedback = profile.delay.get("feedback", 0.3)
        mix = profile.delay.get("mix", 0.2)
        if mix > 0:
            left = apply_delay(left, sr, bpm=bpm,
                               feedback=feedback, mix=mix)
            right = apply_delay(right, sr, bpm=bpm,
                                subdivision=0.375,
                                feedback=feedback, mix=mix)

    # 8. Re-normalize to -1 dBFS
    target_peak = 10.0 ** (-1.0 / 20.0)  # ~0.891
    peak = max(np.max(np.abs(left)), np.max(np.abs(right)), 1e-10)
    if peak > 0:
        gain = target_peak / peak
        left *= gain
        right *= gain

    # 9. Soft-limit
    left = _soft_limit_channel(left)
    right = _soft_limit_channel(right)

    # 10. Fade in/out 15ms
    left = fade_in_out(left, sr, fade_ms=15.0)
    right = fade_in_out(right, sr, fade_ms=15.0)

    return np.column_stack([left, right])
