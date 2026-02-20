"""Configuration loading with 3-level precedence: defaults → file → CLI overrides."""

from __future__ import annotations

import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

if sys.version_info >= (3, 11):
    import tomllib
else:
    import tomli as tomllib


DEFAULT_CONFIG: dict[str, Any] = {
    "audio": {
        "sample_rate": 44100,
        "bit_depth": 16,
        "mp3_bitrate": 192,
    },
    "mapping": {
        "min_bpm": 70,
        "max_bpm": 180,
        "min_duration": 2.0,
        "max_duration": 120.0,
        "mix": {
            "lead_melody": 0.85,
            "bass_loops": 0.70,
            "counter_melodies": 0.55,
            "drone_pads": 0.40,
            "percussion": 0.50,
            "motifs": 0.65,
            "effects": 0.30,
        },
        "instruments": {},
    },
}

CONFIG_SEARCH_PATHS = [
    Path("chemsymphony.toml"),
    Path.home() / ".config" / "chemsymphony" / "config.toml",
]


@dataclass
class Config:
    """Merged configuration object."""

    audio_sample_rate: int = 44100
    audio_bit_depth: int = 16
    audio_mp3_bitrate: int = 192

    min_bpm: int = 70
    max_bpm: int = 180
    min_duration: float = 2.0
    max_duration: float = 120.0

    mix: dict[str, float] = field(default_factory=lambda: {
        "lead_melody": 0.85,
        "bass_loops": 0.70,
        "counter_melodies": 0.55,
        "drone_pads": 0.40,
        "percussion": 0.50,
        "motifs": 0.65,
        "effects": 0.30,
    })

    instruments: dict[str, str] = field(default_factory=dict)

    # CLI overrides (None means "auto-derive")
    bpm: int | None = None
    duration: float | None = None
    key: str | None = None
    seed: int | None = None
    post_processing: str | None = None


def _deep_merge(base: dict, override: dict) -> dict:
    """Recursively merge *override* into *base*."""
    merged = dict(base)
    for k, v in override.items():
        if k in merged and isinstance(merged[k], dict) and isinstance(v, dict):
            merged[k] = _deep_merge(merged[k], v)
        else:
            merged[k] = v
    return merged


def _find_config_file() -> Path | None:
    for p in CONFIG_SEARCH_PATHS:
        if p.is_file():
            return p
    return None


def load_config(
    *,
    config_path: str | Path | None = None,
    bpm: int | None = None,
    duration: float | None = None,
    key: str | None = None,
    seed: int | None = None,
    post_processing: str | None = None,
) -> Config:
    """Load configuration with 3-level precedence.

    1. Built-in defaults
    2. TOML config file (auto-discovered or explicit)
    3. CLI keyword overrides
    """
    merged = dict(DEFAULT_CONFIG)

    # Level 2: TOML file
    file_path = Path(config_path) if config_path else _find_config_file()
    if file_path and file_path.is_file():
        with open(file_path, "rb") as f:
            file_data = tomllib.load(f)
        merged = _deep_merge(merged, file_data)

    audio = merged.get("audio", {})
    mapping = merged.get("mapping", {})

    cfg = Config(
        audio_sample_rate=audio.get("sample_rate", 44100),
        audio_bit_depth=audio.get("bit_depth", 16),
        audio_mp3_bitrate=audio.get("mp3_bitrate", 192),
        min_bpm=mapping.get("min_bpm", 70),
        max_bpm=mapping.get("max_bpm", 180),
        min_duration=mapping.get("min_duration", 2.0),
        max_duration=mapping.get("max_duration", 120.0),
        mix=mapping.get("mix", DEFAULT_CONFIG["mapping"]["mix"]),
        instruments=mapping.get("instruments", {}),
        bpm=bpm,
        duration=duration,
        key=key,
        seed=seed,
        post_processing=post_processing,
    )
    return cfg
