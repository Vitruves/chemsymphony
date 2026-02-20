"""Techno genre profile: driving, repetitive, bass-heavy electronic music."""

from chemsymphony.genres import GenreProfile

GENRE = GenreProfile(
    name="Techno",
    bpm_range=(125, 145),
    mix_adjustments={
        "bass_loops": 1.3,
        "percussion": 1.4,
        "drone_pads": 0.5,
        "lead_melody": 0.8,
    },
    eq={
        "low_boost": 0.5,
        "mid_boost": -0.3,
        "brightness": 0.1,
    },
    compression={
        "threshold_db": -15.0,
        "ratio": 6.0,
    },
    reverb_depth=0.15,
    saturation=0.2,
    highpass_hz=30.0,
)
