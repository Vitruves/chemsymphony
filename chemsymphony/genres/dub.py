"""Dub genre profile: deep bass, spacious delay, heavy reverb."""

from chemsymphony.genres import GenreProfile

GENRE = GenreProfile(
    name="Dub",
    bpm_range=(60, 90),
    mix_adjustments={
        "bass_loops": 1.5,
        "percussion": 1.0,
        "lead_melody": 0.7,
        "drone_pads": 1.2,
        "effects": 1.4,
        "counter_melodies": 0.6,
    },
    eq={
        "low_boost": 0.7,
        "mid_boost": -0.2,
        "brightness": -0.3,
    },
    compression={
        "threshold_db": -10.0,
        "ratio": 3.0,
    },
    reverb_depth=0.6,
    saturation=0.1,
    highpass_hz=25.0,
    delay={"feedback": 0.5, "mix": 0.35},
    swing_override=0.1,
)
