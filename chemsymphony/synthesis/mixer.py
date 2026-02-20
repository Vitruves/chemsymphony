"""Layer mixing and mastering with compression and proper gain staging."""

from __future__ import annotations

import numpy as np
from scipy import signal as scipy_signal

from chemsymphony.config import Config
from chemsymphony.features import MolecularFeatures
from chemsymphony.synthesis.effects import (
    apply_reverb, apply_eq, apply_chorus, apply_delay, fade_in_out,
)


# ---------------------------------------------------------------------------
# Dynamic range compression
# ---------------------------------------------------------------------------

def _compress(audio: np.ndarray, sr: int, threshold_db: float = -12.0,
              ratio: float = 4.0, attack_ms: float = 5.0,
              release_ms: float = 50.0) -> np.ndarray:
    """RMS-based compressor on a mono signal (vectorized envelope)."""
    if len(audio) == 0:
        return audio

    threshold = 10.0 ** (threshold_db / 20.0)
    attack_coeff = np.exp(-1.0 / (attack_ms * sr / 1000.0))
    release_coeff = np.exp(-1.0 / (release_ms * sr / 1000.0))

    # RMS envelope detection with smoothing
    rms_window = int(sr * 0.01)  # 10ms window
    if rms_window < 1:
        rms_window = 1

    # Compute RMS envelope
    squared = audio ** 2
    pad = np.zeros(rms_window)
    padded = np.concatenate([pad, squared])
    cumsum = np.cumsum(padded)
    rms = np.sqrt((cumsum[rms_window:] - cumsum[:-rms_window]) / rms_window)
    rms = rms[:len(audio)]

    # Compute target gain for all samples at once
    target_gain = np.ones_like(audio)
    above = rms > threshold
    if np.any(above):
        compressed_level = threshold + (rms[above] - threshold) / ratio
        target_gain[above] = compressed_level / np.maximum(rms[above], 1e-10)

    # Smooth with asymmetric attack/release envelope follower
    gain = np.empty_like(audio)
    envelope = 1.0
    for i in range(len(audio)):
        tg = target_gain[i]
        if tg < envelope:
            envelope = attack_coeff * envelope + (1.0 - attack_coeff) * tg
        else:
            envelope = release_coeff * envelope + (1.0 - release_coeff) * tg
        gain[i] = envelope

    return audio * gain


def _compress_stereo(left: np.ndarray, right: np.ndarray, sr: int,
                     threshold_db: float = -12.0,
                     ratio: float = 4.0) -> tuple[np.ndarray, np.ndarray]:
    """Apply linked stereo compression (same gain to both channels)."""
    # Use the louder channel for detection
    combined = np.maximum(np.abs(left), np.abs(right))
    # Compress the combined signal to get the gain curve
    compressed = _compress(combined, sr, threshold_db, ratio)
    gain = np.ones_like(combined)
    mask = combined > 1e-10
    gain[mask] = compressed[mask] / combined[mask]

    return left * gain, right * gain


# ---------------------------------------------------------------------------
# Soft limiter
# ---------------------------------------------------------------------------

def _soft_limit(audio: np.ndarray, headroom_db: float = -1.0) -> np.ndarray:
    """Soft limiter using tanh with headroom and 2x oversampling."""
    ceiling = 10.0 ** (headroom_db / 20.0)  # ~0.89 for -1dB

    # 2x upsample
    n = len(audio)
    upsampled = np.zeros(n * 2, dtype=audio.dtype)
    upsampled[::2] = audio
    upsampled[1::2] = audio
    # Smooth interpolation via simple averaging
    upsampled[1:-1:2] = (audio[:-1] + audio[1:]) * 0.5

    # Apply nonlinearity at higher sample rate
    limited = np.tanh(upsampled / ceiling) * ceiling

    # Anti-alias lowpass before downsampling
    sos = scipy_signal.butter(4, 0.5, btype="low", output="sos")
    limited = scipy_signal.sosfilt(sos, limited)

    # 2x downsample
    return limited[::2]


# ---------------------------------------------------------------------------
# Main mixer
# ---------------------------------------------------------------------------

def mix_layers(
    raw_layers: list[tuple[str, np.ndarray, float]],
    feat: MolecularFeatures,
    cfg: Config,
) -> np.ndarray:
    """Mix all synthesized stereo layers into a stereo master.

    Expects each layer's audio to be (N, 2) stereo.
    Returns a 2D numpy array of shape (n_samples, 2) normalized to [-1, 1].

    Gain staging:
    1. Sum layers with volumes → raw mix
    2. Bus compression
    3. Global effects (reverb, EQ, chorus)
    4. Normalize to -1 dBFS
    5. Soft limiter
    6. Fade in/out (15ms)
    7. Trim silence
    """
    sr = cfg.audio_sample_rate
    ap = feat.audio_parameters
    duration_sec = ap["duration"]
    total_samples = int(sr * (duration_sec + 1))

    left = np.zeros(total_samples, dtype=np.float64)
    right = np.zeros(total_samples, dtype=np.float64)

    for name, audio, volume in raw_layers:
        # Audio is stereo (N, 2)
        if audio.ndim == 2 and audio.shape[1] == 2:
            l_ch = audio[:, 0]
            r_ch = audio[:, 1]
        else:
            # Fallback for mono: apply center pan
            l_ch = audio * 0.707
            r_ch = audio * 0.707

        # Trim or pad to total_samples
        n = len(l_ch)
        if n > total_samples:
            l_ch = l_ch[:total_samples]
            r_ch = r_ch[:total_samples]
        elif n < total_samples:
            l_pad = np.zeros(total_samples)
            r_pad = np.zeros(total_samples)
            l_pad[:n] = l_ch
            r_pad[:n] = r_ch
            l_ch = l_pad
            r_ch = r_pad

        left += l_ch * volume
        right += r_ch * volume

    # Step 1.5: DC offset removal (5Hz highpass, Butterworth order 1)
    dc_cutoff = 5.0 / (sr / 2)
    if 0 < dc_cutoff < 1:
        dc_sos = scipy_signal.butter(1, dc_cutoff, btype="high", output="sos")
        left = scipy_signal.sosfilt(dc_sos, left)
        right = scipy_signal.sosfilt(dc_sos, right)

    # Step 2: Bus compression — adaptive threshold and ratio from molecular complexity
    # Higher bertz_ct (more complex) → more aggressive compression (lower threshold)
    # Higher aromatic_fraction → tighter ratio (more controlled dynamics)
    bertz_ct = feat.bertz_ct
    aromatic_frac = feat.aromatic_fraction
    comp_threshold = -9.0 - 6.0 * min(bertz_ct / 1500.0, 1.0)   # -9 to -15 dB
    comp_ratio = 2.0 + 3.0 * aromatic_frac                        # 2:1 to 5:1
    left, right = _compress_stereo(left, right, sr,
                                   threshold_db=comp_threshold,
                                   ratio=comp_ratio)

    # Step 3: Global effects
    # Reverb — use reverb_wetness from audio params (falls back to aromatic count)
    reverb_depth = ap.get("reverb_wetness", min(1.0, feat.aromatic_atom_count / 20.0))
    reverb_diffusion = ap.get("reverb_diffusion", 0.0)
    # Diffusion increases damping (more diffuse = smoother tail)
    reverb_damping = 0.3 + 0.4 * reverb_diffusion  # 0.3–0.7
    if reverb_depth > 0.05:
        left = apply_reverb(left, sr, depth=reverb_depth, damping=reverb_damping)
        right = apply_reverb(right, sr, depth=reverb_depth, damping=reverb_damping)

    # EQ from net charge
    brightness = 0.0
    if feat.net_charge > 0:
        brightness = min(0.5, feat.net_charge * 0.2)
    elif feat.net_charge < 0:
        brightness = max(-0.5, feat.net_charge * 0.2)

    # Also use filter_warmth for low-end boost
    filter_warmth = ap.get("filter_warmth", 0.0)
    low_boost = filter_warmth * 0.3

    if abs(brightness) > 0.05 or low_boost > 0.05:
        left = apply_eq(left, sr, brightness=brightness, low_boost=low_boost)
        right = apply_eq(right, sr, brightness=brightness, low_boost=low_boost)

    # Chorus for organic timbres
    timbre_organic = ap.get("timbre_organic", 0.0)
    if timbre_organic > 0.4:
        chorus_mix = (timbre_organic - 0.4) * 0.5
        left = apply_chorus(left, sr, rate=0.8, depth=0.002, voices=2, mix=chorus_mix)
        right = apply_chorus(right, sr, rate=0.8, depth=0.002, voices=2, mix=chorus_mix)

    # Echo density: terminal atoms → delay reflections
    echo_density = ap.get("echo_density", 0.0)
    if echo_density > 0.1:
        delay_mix = echo_density * 0.2  # 0–0.2
        delay_feedback = 0.2 + echo_density * 0.2  # 0.2–0.4
        left = apply_delay(left, sr, bpm=ap["bpm"], subdivision=0.25,
                           feedback=delay_feedback, mix=delay_mix)
        right = apply_delay(right, sr, bpm=ap["bpm"], subdivision=0.375,
                            feedback=delay_feedback, mix=delay_mix)

    # Filter sweep: EN variance → slow filter modulation via EQ mid boost
    filter_sweep = ap.get("filter_sweep", 0.0)
    if filter_sweep > 0.1:
        mid_boost = filter_sweep * 0.3  # 0–0.3 mid emphasis
        left = apply_eq(left, sr, mid_boost=mid_boost)
        right = apply_eq(right, sr, mid_boost=mid_boost)

    # Step 4: Normalize to -1 dBFS (0.89 peak)
    target_peak = 10.0 ** (-1.0 / 20.0)  # ~0.891
    peak = max(np.max(np.abs(left)), np.max(np.abs(right)), 1e-10)
    if peak > 0:
        gain = target_peak / peak
        left *= gain
        right *= gain

    # Step 5: Soft limiter
    left = _soft_limit(left)
    right = _soft_limit(right)

    # Step 6: Fade in/out (15ms)
    left = fade_in_out(left, sr, fade_ms=15.0)
    right = fade_in_out(right, sr, fade_ms=15.0)

    # Step 7: Trim silence from end
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
