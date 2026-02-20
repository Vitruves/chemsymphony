"""Ambient genre profile: slow, spacious, texture-rich atmospheric music."""

from chemsymphony.genres import GenreProfile

GENRE = GenreProfile(
    name="Ambient",
    bpm_range=(60, 90),
    mix_adjustments={
        "drone_pads": 1.8,
        "effects": 1.5,
        "percussion": 0.2,
        "bass_loops": 0.6,
        "lead_melody": 0.7,
        "counter_melodies": 1.2,
    },
    eq={
        "low_boost": 0.2,
        "mid_boost": -0.1,
        "brightness": -0.2,
    },
    compression={
        "threshold_db": -8.0,
        "ratio": 2.0,
    },
    reverb_depth=0.7,
    saturation=0.0,
    lowpass_hz=12000.0,
    delay={"feedback": 0.4, "mix": 0.3},
)
