"""House genre profile: warm, groovy, four-on-the-floor electronic music."""

from chemsymphony.genres import GenreProfile

GENRE = GenreProfile(
    name="House",
    bpm_range=(118, 130),
    mix_adjustments={
        "bass_loops": 1.2,
        "percussion": 1.1,
        "drone_pads": 0.8,
        "counter_melodies": 1.1,
    },
    eq={
        "low_boost": 0.4,
        "mid_boost": 0.1,
        "brightness": -0.1,
    },
    compression={
        "threshold_db": -12.0,
        "ratio": 4.0,
    },
    reverb_depth=0.3,
    saturation=0.1,
    highpass_hz=25.0,
    swing_override=0.15,
)
