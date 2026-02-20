"""Drum & bass genre profile: fast breakbeats, deep sub-bass, high energy."""

from chemsymphony.genres import GenreProfile

GENRE = GenreProfile(
    name="Drum & Bass",
    bpm_range=(170, 180),
    mix_adjustments={
        "percussion": 1.5,
        "bass_loops": 1.4,
        "lead_melody": 0.9,
        "drone_pads": 0.4,
        "counter_melodies": 0.7,
    },
    eq={
        "low_boost": 0.6,
        "mid_boost": -0.2,
        "brightness": 0.3,
    },
    compression={
        "threshold_db": -16.0,
        "ratio": 6.0,
    },
    reverb_depth=0.15,
    saturation=0.25,
    highpass_hz=30.0,
)
