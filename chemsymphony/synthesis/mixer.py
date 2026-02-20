"""Layer mixing and mastering."""

from __future__ import annotations

import numpy as np

from chemsymphony.config import Config
from chemsymphony.features import MolecularFeatures
from chemsymphony.synthesis.effects import apply_reverb, apply_eq, fade_in_out


def mix_layers(
    raw_layers: list[tuple[str, np.ndarray, float]],
    feat: MolecularFeatures,
    cfg: Config,
) -> np.ndarray:
    """Mix all synthesized layers into a stereo master.

    Returns a 2D numpy array of shape (n_samples, 2) normalized to [-1, 1].
    """
    sr = cfg.audio_sample_rate
    ap = feat.audio_parameters
    duration_sec = ap["duration"]
    total_samples = int(sr * (duration_sec + 1))  # +1s for tails

    # Stereo output
    left = np.zeros(total_samples, dtype=np.float64)
    right = np.zeros(total_samples, dtype=np.float64)

    for name, audio, volume in raw_layers:
        # Ensure same length
        if len(audio) > total_samples:
            audio = audio[:total_samples]
        elif len(audio) < total_samples:
            padded = np.zeros(total_samples)
            padded[:len(audio)] = audio
            audio = padded

        # Apply volume
        audio = audio * volume

        # Simple pan: center for now (expression could refine this)
        left += audio * 0.707
        right += audio * 0.707

    # Apply global effects
    from chemsymphony.mapping.expression import ExpressionParams

    # Apply reverb based on aromatic content
    reverb_depth = min(1.0, feat.aromatic_atom_count / 20.0)
    if reverb_depth > 0.05:
        left = apply_reverb(left, sr, depth=reverb_depth)
        right = apply_reverb(right, sr, depth=reverb_depth)

    # Apply EQ based on net charge
    brightness = 0.0
    if feat.net_charge > 0:
        brightness = min(0.5, feat.net_charge * 0.2)
    elif feat.net_charge < 0:
        brightness = max(-0.5, feat.net_charge * 0.2)
    if abs(brightness) > 0.05:
        left = apply_eq(left, sr, brightness)
        right = apply_eq(right, sr, brightness)

    # Normalize
    peak = max(np.max(np.abs(left)), np.max(np.abs(right)), 1e-10)
    if peak > 1.0:
        left /= peak
        right /= peak

    # Soft limit
    left = np.tanh(left)
    right = np.tanh(right)

    # Fade in/out to remove clicks
    left = fade_in_out(left, sr)
    right = fade_in_out(right, sr)

    # Trim silence from end
    threshold = 0.001
    combined = np.abs(left) + np.abs(right)
    nonzero = np.where(combined > threshold)[0]
    if len(nonzero) > 0:
        end = min(nonzero[-1] + sr, total_samples)  # Keep 1s tail
    else:
        end = int(sr * duration_sec)
    left = left[:end]
    right = right[:end]

    return np.column_stack([left, right])
