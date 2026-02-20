"""Trance genre profile: euphoric, layered pads, driving arpeggios."""

from chemsymphony.genres import GenreProfile

GENRE = GenreProfile(
    name="Trance",
    bpm_range=(138, 150),
    mix_adjustments={
        "drone_pads": 1.4,
        "lead_melody": 1.2,
        "counter_melodies": 1.2,
        "bass_loops": 1.0,
        "percussion": 0.9,
    },
    eq={
        "low_boost": 0.3,
        "mid_boost": 0.2,
        "brightness": 0.3,
    },
    compression={
        "threshold_db": -14.0,
        "ratio": 5.0,
    },
    reverb_depth=0.45,
    saturation=0.1,
    highpass_hz=25.0,
    delay={"feedback": 0.35, "mix": 0.2},
)
