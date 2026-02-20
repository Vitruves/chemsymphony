"""Audio effects: reverb, panning, EQ, distortion, vibrato."""

from __future__ import annotations

import numpy as np
from scipy import signal as scipy_signal


def apply_reverb(audio: np.ndarray, sr: int, depth: float = 0.3,
                 predelay_ms: float = 20.0) -> np.ndarray:
    """Simple Schroeder reverb using comb + allpass filters."""
    if depth <= 0:
        return audio

    depth = min(1.0, depth)
    output = audio.copy()

    predelay_samples = int(predelay_ms * sr / 1000)
    if predelay_samples > 0 and predelay_samples < len(audio):
        delayed = np.zeros_like(audio)
        delayed[predelay_samples:] = audio[:-predelay_samples]
    else:
        delayed = audio

    # Comb filters at different delays
    comb_delays_ms = [29.7, 37.1, 41.1, 43.7]
    comb_gains = [0.805, 0.827, 0.783, 0.764]

    for delay_ms, gain in zip(comb_delays_ms, comb_gains):
        delay_samples = int(delay_ms * sr / 1000)
        if delay_samples >= len(delayed):
            continue
        comb = np.zeros_like(delayed)
        for i in range(delay_samples, len(delayed)):
            comb[i] = delayed[i] + gain * comb[i - delay_samples]
        output += comb * depth * 0.25

    return output


def apply_stereo_pan(left: np.ndarray, right: np.ndarray, pan: float) -> tuple[np.ndarray, np.ndarray]:
    """Apply stereo panning. pan: -1.0 (left) to 1.0 (right)."""
    pan = max(-1.0, min(1.0, pan))
    l_gain = np.cos((pan + 1) * np.pi / 4)
    r_gain = np.sin((pan + 1) * np.pi / 4)
    return left * l_gain, right * r_gain


def apply_eq(audio: np.ndarray, sr: int, brightness: float = 0.0) -> np.ndarray:
    """3-band EQ. brightness: -1.0 (dark) to 1.0 (bright)."""
    if abs(brightness) < 0.05:
        return audio

    if brightness > 0:
        # Boost highs
        b, a = scipy_signal.butter(2, 4000 / (sr / 2), btype="high")
        high = scipy_signal.lfilter(b, a, audio)
        return audio + high * brightness * 0.3
    else:
        # Boost lows
        b, a = scipy_signal.butter(2, 300 / (sr / 2), btype="low")
        low = scipy_signal.lfilter(b, a, audio)
        return audio + low * abs(brightness) * 0.3


def apply_soft_distortion(audio: np.ndarray, amount: float = 0.5) -> np.ndarray:
    """Soft clipping via tanh."""
    if amount <= 0:
        return audio
    gain = 1.0 + amount * 4.0
    return np.tanh(audio * gain)


def apply_vibrato(audio: np.ndarray, sr: int, rate: float = 5.0,
                  depth: float = 0.002) -> np.ndarray:
    """Apply vibrato (pitch modulation) using variable delay."""
    if depth <= 0:
        return audio
    n = len(audio)
    t = np.arange(n) / sr
    mod = depth * sr * np.sin(2 * np.pi * rate * t)
    indices = np.arange(n) + mod
    indices = np.clip(indices, 0, n - 1).astype(int)
    return audio[indices]


def fade_in_out(audio: np.ndarray, sr: int, fade_ms: float = 2.0) -> np.ndarray:
    """Apply fade in/out to prevent clicks."""
    fade_samples = int(fade_ms * sr / 1000)
    fade_samples = min(fade_samples, len(audio) // 2)
    if fade_samples <= 0:
        return audio
    out = audio.copy()
    out[:fade_samples] *= np.linspace(0, 1, fade_samples)
    out[-fade_samples:] *= np.linspace(1, 0, fade_samples)
    return out
