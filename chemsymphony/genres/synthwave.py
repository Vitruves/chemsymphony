"""Synthwave genre profile: retro 80s synths, warm, nostalgic."""

from chemsymphony.genres import GenreProfile

GENRE = GenreProfile(
    name="Synthwave",
    bpm_range=(80, 118),
    mix_adjustments={
        "lead_melody": 1.2,
        "drone_pads": 1.3,
        "bass_loops": 1.1,
        "percussion": 0.8,
        "counter_melodies": 1.1,
    },
    eq={
        "low_boost": 0.4,
        "mid_boost": 0.1,
        "brightness": 0.2,
    },
    compression={
        "threshold_db": -12.0,
        "ratio": 4.0,
    },
    reverb_depth=0.4,
    saturation=0.15,
    delay={"feedback": 0.3, "mix": 0.2},
)
