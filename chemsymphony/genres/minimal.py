"""Minimal genre profile: sparse, clean, subtle micro-variations."""

from chemsymphony.genres import GenreProfile

GENRE = GenreProfile(
    name="Minimal",
    bpm_range=(120, 135),
    mix_adjustments={
        "percussion": 1.2,
        "bass_loops": 1.1,
        "lead_melody": 0.6,
        "drone_pads": 0.4,
        "counter_melodies": 0.5,
        "motifs": 0.7,
        "effects": 0.8,
    },
    eq={
        "low_boost": 0.2,
        "mid_boost": -0.1,
        "brightness": 0.1,
    },
    compression={
        "threshold_db": -10.0,
        "ratio": 3.0,
    },
    reverb_depth=0.25,
    saturation=0.05,
    highpass_hz=30.0,
)
