"""Industrial genre profile: harsh, metallic, aggressive electronic music."""

from chemsymphony.genres import GenreProfile

GENRE = GenreProfile(
    name="Industrial",
    bpm_range=(120, 140),
    mix_adjustments={
        "percussion": 1.5,
        "bass_loops": 1.2,
        "lead_melody": 1.1,
        "drone_pads": 0.6,
        "effects": 1.3,
    },
    eq={
        "low_boost": 0.3,
        "mid_boost": 0.4,
        "brightness": 0.5,
    },
    compression={
        "threshold_db": -18.0,
        "ratio": 8.0,
    },
    reverb_depth=0.2,
    saturation=0.55,
    highpass_hz=35.0,
)
