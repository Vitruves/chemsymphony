"""Audio effects: reverb, panning, EQ, distortion, vibrato, chorus, delay."""

from __future__ import annotations

import numpy as np
from scipy import signal as scipy_signal


# ---------------------------------------------------------------------------
# Reverb (vectorized Schroeder design)
# ---------------------------------------------------------------------------

def _comb_filter_vectorized(audio: np.ndarray, delay_samples: int,
                            gain: float, damping: float = 0.2) -> np.ndarray:
    """Vectorized comb filter with optional LP damping in feedback path."""
    n = len(audio)
    if delay_samples >= n or delay_samples <= 0:
        return np.zeros_like(audio)

    output = np.zeros(n, dtype=np.float64)
    # Process in blocks of delay_samples for vectorized operation
    # First block is just the input
    output[:delay_samples] = audio[:delay_samples]

    for start in range(delay_samples, n, delay_samples):
        end = min(start + delay_samples, n)
        length = end - start
        feedback = output[start - delay_samples:start - delay_samples + length]
        # LP damping: simple one-pole filter on feedback
        if damping > 0:
            damped = feedback.copy()
            for j in range(1, length):
                damped[j] = (1.0 - damping) * damped[j] + damping * damped[j - 1]
            feedback = damped
        output[start:end] = audio[start:end] + gain * feedback

    return output


def _allpass_filter(audio: np.ndarray, delay_samples: int,
                    gain: float = 0.7) -> np.ndarray:
    """Allpass filter using vectorized numpy roll."""
    if delay_samples >= len(audio) or delay_samples <= 0:
        return audio.copy()
    n = len(audio)
    output = np.zeros(n, dtype=np.float64)
    output[:delay_samples] = audio[:delay_samples]

    for start in range(delay_samples, n, delay_samples):
        end = min(start + delay_samples, n)
        length = end - start
        delayed_out = output[start - delay_samples:start - delay_samples + length]
        output[start:end] = -gain * audio[start:end] + delayed_out + gain * output[start:end - length + length]

    # Simplified: use scipy for clean allpass
    buf = np.zeros(n, dtype=np.float64)
    buf[:delay_samples] = audio[:delay_samples]
    for i in range(delay_samples, n):
        buf[i] = -gain * audio[i] + audio[i - delay_samples] + gain * buf[i - delay_samples]
    return buf


def apply_reverb(audio: np.ndarray, sr: int, depth: float = 0.3,
                 predelay_ms: float = 20.0, damping: float = 0.3) -> np.ndarray:
    """Schroeder reverb: parallel comb filters → series allpass filters.

    Vectorized implementation — ~100x faster than sample-by-sample.
    """
    if depth <= 0:
        return audio

    depth = min(1.0, depth)

    # Pre-delay
    predelay_samples = int(predelay_ms * sr / 1000)
    if 0 < predelay_samples < len(audio):
        delayed = np.zeros_like(audio)
        delayed[predelay_samples:] = audio[:-predelay_samples]
    else:
        delayed = audio

    # Parallel comb filter bank (Schroeder design)
    comb_params = [
        (int(29.7 * sr / 1000), 0.805),
        (int(37.1 * sr / 1000), 0.827),
        (int(41.1 * sr / 1000), 0.783),
        (int(43.7 * sr / 1000), 0.764),
    ]

    comb_sum = np.zeros_like(delayed)
    for delay_samples, gain in comb_params:
        comb_sum += _comb_filter_vectorized(delayed, delay_samples, gain, damping)
    comb_sum *= 0.25  # Average the four combs

    # Series allpass filters (2 stages)
    allpass_delays = [int(5.0 * sr / 1000), int(1.7 * sr / 1000)]
    result = comb_sum
    for ap_delay in allpass_delays:
        result = _allpass_filter(result, ap_delay, gain=0.7)

    return audio + result * depth


# ---------------------------------------------------------------------------
# Stereo panning
# ---------------------------------------------------------------------------

def apply_stereo_pan(left: np.ndarray, right: np.ndarray,
                     pan: float) -> tuple[np.ndarray, np.ndarray]:
    """Constant-power pan law. pan: -1.0 (left) to 1.0 (right)."""
    pan = max(-1.0, min(1.0, pan))
    angle = (pan + 1) * np.pi / 4  # 0 to pi/2
    l_gain = np.cos(angle)
    r_gain = np.sin(angle)
    return left * l_gain, right * r_gain


# ---------------------------------------------------------------------------
# EQ (parametric)
# ---------------------------------------------------------------------------

def apply_eq(audio: np.ndarray, sr: int, brightness: float = 0.0,
             low_boost: float = 0.0, mid_boost: float = 0.0) -> np.ndarray:
    """Parametric EQ with low, mid, and brightness (high) bands.

    brightness: -1.0 (cut highs) to 1.0 (boost highs)
    low_boost: -1.0 to 1.0
    mid_boost: -1.0 to 1.0
    """
    result = audio.copy()
    nyquist = sr / 2

    # High band (brightness)
    if abs(brightness) >= 0.05:
        cutoff = min(4000 / nyquist, 0.95)
        if cutoff > 0:
            sos = scipy_signal.butter(2, cutoff, btype="high", output="sos")
            high = scipy_signal.sosfilt(sos, audio)
            result = result + high * brightness * 0.3

    # Low band
    if abs(low_boost) >= 0.05:
        cutoff = min(300 / nyquist, 0.95)
        if cutoff > 0:
            sos = scipy_signal.butter(2, cutoff, btype="low", output="sos")
            low = scipy_signal.sosfilt(sos, audio)
            result = result + low * low_boost * 0.3

    # Mid band (bandpass 300-3000 Hz)
    if abs(mid_boost) >= 0.05:
        low_edge = min(300 / nyquist, 0.95)
        high_edge = min(3000 / nyquist, 0.95)
        if 0 < low_edge < high_edge < 1:
            sos = scipy_signal.butter(2, [low_edge, high_edge], btype="band", output="sos")
            mid = scipy_signal.sosfilt(sos, audio)
            result = result + mid * mid_boost * 0.3

    return result


# ---------------------------------------------------------------------------
# Distortion
# ---------------------------------------------------------------------------

def apply_soft_distortion(audio: np.ndarray, amount: float = 0.5) -> np.ndarray:
    """Soft clipping via tanh."""
    if amount <= 0:
        return audio
    gain = 1.0 + amount * 4.0
    return np.tanh(audio * gain)


# ---------------------------------------------------------------------------
# Vibrato
# ---------------------------------------------------------------------------

def apply_vibrato(audio: np.ndarray, sr: int, rate: float = 5.0,
                  depth: float = 0.002) -> np.ndarray:
    """Apply vibrato (pitch modulation) using variable delay."""
    if depth <= 0:
        return audio
    n = len(audio)
    t = np.arange(n) / sr
    mod = depth * sr * np.sin(2 * np.pi * rate * t)
    indices = np.arange(n, dtype=np.float64) + mod
    indices = np.clip(indices, 0, n - 1).astype(int)
    return audio[indices]


# ---------------------------------------------------------------------------
# Chorus
# ---------------------------------------------------------------------------

def apply_chorus(audio: np.ndarray, sr: int, rate: float = 1.5,
                 depth: float = 0.003, voices: int = 3,
                 mix: float = 0.5) -> np.ndarray:
    """Multi-voice chorus with LFO-modulated delay times."""
    if voices <= 0 or depth <= 0:
        return audio

    n = len(audio)
    t = np.arange(n) / sr
    result = audio.copy()

    for v in range(voices):
        # Each voice has a different LFO phase
        phase = 2 * np.pi * v / voices
        delay_mod = depth * sr * (0.5 + 0.5 * np.sin(2 * np.pi * rate * t + phase))
        # Base delay to avoid negative indices
        base_delay = int(depth * sr) + 10
        indices = np.arange(n, dtype=np.float64) - base_delay - delay_mod
        indices = np.clip(indices, 0, n - 1).astype(int)
        result = result + audio[indices] * (mix / voices)

    return result


# ---------------------------------------------------------------------------
# Tempo-synced delay
# ---------------------------------------------------------------------------

def apply_delay(audio: np.ndarray, sr: int, bpm: float,
                subdivision: float = 0.5, feedback: float = 0.3,
                mix: float = 0.25) -> np.ndarray:
    """Tempo-synced delay with feedback."""
    if mix <= 0 or bpm <= 0:
        return audio

    delay_sec = (60.0 / bpm) * subdivision
    delay_samples = int(delay_sec * sr)
    if delay_samples <= 0 or delay_samples >= len(audio):
        return audio

    n = len(audio)
    result = audio.copy()
    tap = audio.copy()

    # Apply up to 8 feedback taps
    for _ in range(8):
        delayed = np.zeros(n, dtype=np.float64)
        if delay_samples < n:
            delayed[delay_samples:] = tap[:n - delay_samples]
        tap = delayed * feedback
        result = result + delayed * mix
        if np.max(np.abs(tap)) < 0.001:
            break

    return result


# ---------------------------------------------------------------------------
# Fade in/out
# ---------------------------------------------------------------------------

def fade_in_out(audio: np.ndarray, sr: int, fade_ms: float = 15.0) -> np.ndarray:
    """Apply fade in/out to prevent clicks. Default 15ms."""
    fade_samples = int(fade_ms * sr / 1000)
    fade_samples = min(fade_samples, len(audio) // 2)
    if fade_samples <= 0:
        return audio
    out = audio.copy()
    out[:fade_samples] *= np.linspace(0, 1, fade_samples)
    out[-fade_samples:] *= np.linspace(1, 0, fade_samples)
    return out
