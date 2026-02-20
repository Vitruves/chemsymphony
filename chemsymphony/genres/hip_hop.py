"""Hip hop genre profile: boom-bap, heavy bass, swung grooves."""

from chemsymphony.genres import GenreProfile

GENRE = GenreProfile(
    name="Hip Hop",
    bpm_range=(80, 100),
    mix_adjustments={
        "bass_loops": 1.4,
        "percussion": 1.3,
        "lead_melody": 0.9,
        "drone_pads": 0.5,
        "motifs": 1.1,
    },
    eq={
        "low_boost": 0.6,
        "mid_boost": -0.1,
        "brightness": -0.15,
    },
    compression={
        "threshold_db": -14.0,
        "ratio": 5.0,
    },
    reverb_depth=0.2,
    saturation=0.15,
    highpass_hz=30.0,
    swing_override=0.2,
)
