"""Acid techno genre profile: aggressive, squelchy, high-energy electronic music."""

from chemsymphony.genres import GenreProfile

GENRE = GenreProfile(
    name="Acid Techno",
    bpm_range=(130, 150),
    mix_adjustments={
        "bass_loops": 1.4,
        "percussion": 1.3,
        "lead_melody": 1.2,
        "drone_pads": 0.3,
    },
    eq={
        "low_boost": 0.4,
        "mid_boost": 0.2,
        "brightness": 0.6,
    },
    compression={
        "threshold_db": -18.0,
        "ratio": 8.0,
    },
    reverb_depth=0.1,
    saturation=0.6,
    highpass_hz=35.0,
)
