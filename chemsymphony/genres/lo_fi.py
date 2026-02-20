"""Lo-fi genre profile: warm, filtered, vinyl-crackle aesthetic."""

from chemsymphony.genres import GenreProfile

GENRE = GenreProfile(
    name="Lo-Fi",
    bpm_range=(70, 90),
    mix_adjustments={
        "drone_pads": 1.3,
        "lead_melody": 0.8,
        "bass_loops": 1.0,
        "percussion": 0.7,
        "counter_melodies": 1.1,
        "effects": 1.2,
    },
    eq={
        "low_boost": 0.3,
        "mid_boost": 0.15,
        "brightness": -0.4,
    },
    compression={
        "threshold_db": -10.0,
        "ratio": 3.0,
    },
    reverb_depth=0.3,
    saturation=0.15,
    lowpass_hz=10000.0,
    swing_override=0.12,
)
