"""Gabber genre profile: extreme tempo, heavy distortion, aggressive electronic music."""

from chemsymphony.genres import GenreProfile

GENRE = GenreProfile(
    name="Gabber",
    bpm_range=(160, 200),
    mix_adjustments={
        "percussion": 1.6,
        "bass_loops": 1.4,
        "lead_melody": 1.0,
        "drone_pads": 0.2,
        "counter_melodies": 0.5,
        "effects": 0.4,
    },
    eq={
        "low_boost": 0.6,
        "mid_boost": 0.3,
        "brightness": 0.4,
    },
    compression={
        "threshold_db": -20.0,
        "ratio": 10.0,
    },
    reverb_depth=0.05,
    saturation=0.8,
    highpass_hz=40.0,
)
