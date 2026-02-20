"""Jazz genre profile: warm, swinging, harmonically rich acoustic music."""

from chemsymphony.genres import GenreProfile

GENRE = GenreProfile(
    name="Jazz",
    bpm_range=(80, 140),
    scale_preferences=["dorian", "mixolydian", "blues", "bebop_dominant"],
    mix_adjustments={
        "counter_melodies": 1.3,
        "bass_loops": 1.1,
        "lead_melody": 1.1,
        "drone_pads": 0.5,
        "percussion": 0.8,
    },
    eq={
        "low_boost": 0.3,
        "mid_boost": 0.2,
        "brightness": -0.2,
    },
    compression={
        "threshold_db": -8.0,
        "ratio": 2.5,
    },
    reverb_depth=0.35,
    saturation=0.05,
    lowpass_hz=14000.0,
    swing_override=0.25,
)
