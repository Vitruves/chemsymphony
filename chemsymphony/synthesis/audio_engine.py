"""Direct audio synthesis engine using numpy/scipy.

Includes band-limited oscillators (PolyBLEP), resonant filters,
improved instrument models, stereo output, and humanization.
"""

from __future__ import annotations

import numpy as np
from scipy import signal as scipy_signal

from chemsymphony.config import Config
from chemsymphony.features import MolecularFeatures
from chemsymphony.mapping import CompositionLayers
from chemsymphony.mapping.melody import Layer, NoteEvent

# ---------------------------------------------------------------------------
# PolyBLEP correction
# ---------------------------------------------------------------------------

def _polyblep(t: np.ndarray, dt: float) -> np.ndarray:
    """PolyBLEP correction for band-limiting discontinuities.

    *t* is the phase (0–1), *dt* is the phase increment per sample.
    """
    out = np.zeros_like(t)

    # Rising edge (t near 0)
    mask = t < dt
    if np.any(mask):
        x = t[mask] / dt
        out[mask] = 2.0 * x - x * x - 1.0

    # Falling edge (t near 1)
    mask = t > 1.0 - dt
    if np.any(mask):
        x = (t[mask] - 1.0) / dt
        out[mask] = x * x + 2.0 * x + 1.0

    return out


# ---------------------------------------------------------------------------
# Oscillators
# ---------------------------------------------------------------------------

def _sine(freq: float, duration: float, sr: int) -> np.ndarray:
    t = np.arange(int(sr * duration)) / sr
    return np.sin(2 * np.pi * freq * t)


def _square(freq: float, duration: float, sr: int) -> np.ndarray:
    """Band-limited square wave using PolyBLEP."""
    n = int(sr * duration)
    dt = freq / sr
    phase = (np.arange(n) * dt) % 1.0
    # Naive square
    sig = np.where(phase < 0.5, 1.0, -1.0)
    # Apply PolyBLEP at both transitions
    sig += _polyblep(phase, dt)
    sig -= _polyblep((phase + 0.5) % 1.0, dt)
    return sig


def _sawtooth(freq: float, duration: float, sr: int) -> np.ndarray:
    """Band-limited sawtooth wave using PolyBLEP."""
    n = int(sr * duration)
    dt = freq / sr
    phase = (np.arange(n) * dt) % 1.0
    # Naive sawtooth
    sig = 2.0 * phase - 1.0
    # Apply PolyBLEP at discontinuity
    sig -= _polyblep(phase, dt)
    return sig


def _triangle(freq: float, duration: float, sr: int) -> np.ndarray:
    """Triangle wave (inherently band-limited from integrated square)."""
    # Integrate band-limited square for a clean triangle
    sq = _square(freq, duration, sr)
    tri = np.cumsum(sq) / sr
    # Normalize
    peak = np.max(np.abs(tri))
    if peak > 0:
        tri /= peak
    return tri


# ---------------------------------------------------------------------------
# Resonant low-pass filter
# ---------------------------------------------------------------------------

def _lowpass(signal: np.ndarray, cutoff: float, sr: int,
             resonance: float = 0.0, order: int = 2) -> np.ndarray:
    """Apply resonant Butterworth low-pass filter using sosfilt.

    cutoff: frequency in Hz
    resonance: 0.0 (none) to 1.0 (high Q)
    """
    nyquist = sr / 2
    freq = min(cutoff / nyquist, 0.99)
    freq = max(freq, 0.001)

    # Map resonance to a narrower bandwidth (higher order + gain at cutoff)
    actual_order = order + (2 if resonance > 0.5 else 0)
    sos = scipy_signal.butter(actual_order, freq, btype="low", output="sos")
    filtered = scipy_signal.sosfilt(sos, signal)

    # Add resonance peak by mixing in a bandpass at the cutoff
    if resonance > 0.1:
        bw_low = max(freq * 0.8, 0.001)
        bw_high = min(freq * 1.2, 0.999)
        if bw_low < bw_high:
            sos_bp = scipy_signal.butter(2, [bw_low, bw_high], btype="band", output="sos")
            bp = scipy_signal.sosfilt(sos_bp, signal)
            filtered += bp * resonance * 0.5

    return filtered


def _lowpass_stateful(signal: np.ndarray, cutoff: float, sr: int,
                      zi: np.ndarray | None = None,
                      order: int = 2) -> tuple[np.ndarray, np.ndarray]:
    """Stateful Butterworth low-pass filter. Returns (filtered, final_zi).

    Pass *zi* from the previous block to maintain continuity across blocks.
    """
    nyquist = sr / 2
    freq = min(cutoff / nyquist, 0.99)
    freq = max(freq, 0.001)

    sos = scipy_signal.butter(order, freq, btype="low", output="sos")
    if zi is None:
        zi = scipy_signal.sosfilt_zi(sos) * signal[0] if len(signal) > 0 else scipy_signal.sosfilt_zi(sos) * 0.0
    filtered, zo = scipy_signal.sosfilt(sos, signal, zi=zi)
    return filtered, zo


# ---------------------------------------------------------------------------
# ADSR Envelope
# ---------------------------------------------------------------------------

def _adsr(n_samples: int, sr: int, attack: float = 0.01, decay: float = 0.05,
          sustain: float = 0.7, release: float = 0.05) -> np.ndarray:
    """Generate ADSR envelope, clamped to fit within n_samples."""
    if n_samples <= 0:
        return np.zeros(0)

    a = int(attack * sr)
    d = int(decay * sr)
    r = int(release * sr)

    total_adr = a + d + r
    if total_adr > n_samples:
        scale = n_samples / total_adr
        a = int(a * scale)
        d = int(d * scale)
        r = n_samples - a - d

    s = max(0, n_samples - a - d - r)

    env = np.zeros(n_samples)
    pos = 0
    if a > 0:
        env[pos:pos + a] = np.linspace(0, 1, a)
        pos += a
    if d > 0:
        env[pos:pos + d] = sustain + (1.0 - sustain) * np.exp(-5.0 * np.linspace(0, 1, d))
        pos += d
    if s > 0:
        env[pos:pos + s] = sustain
        pos += s
    if r > 0:
        env[pos:pos + r] = sustain * np.exp(-5.0 * np.linspace(0, 1, r))
    return env


# ---------------------------------------------------------------------------
# Instruments
# ---------------------------------------------------------------------------

def _midi_to_freq(midi_note: int) -> float:
    return 440.0 * (2.0 ** ((midi_note - 69) / 12.0))


def _synth_piano(freq: float, duration: float, sr: int,
                 cutoff_scale: float = 1.0) -> np.ndarray:
    """Piano: 8 harmonics with per-harmonic decay, slight chorus detuning."""
    n = int(sr * duration)
    t = np.arange(n) / sr
    sig = np.zeros(n)

    harmonic_amps = [1.0, 0.5, 0.35, 0.25, 0.15, 0.1, 0.06, 0.03]
    harmonic_decays = [2.0, 3.0, 4.0, 5.0, 7.0, 9.0, 11.0, 14.0]

    for k, (amp, decay_rate) in enumerate(zip(harmonic_amps, harmonic_decays), start=1):
        # Slight chorus detuning (±2 cents on odd harmonics)
        detune = 1.0
        if k % 2 == 1 and k > 1:
            detune = 2.0 ** (2.0 / 1200.0)
        sig += amp * np.sin(2 * np.pi * freq * k * detune * t) * np.exp(-decay_rate * 0.3 * t)

    # Velocity-sensitive brightness: higher cutoff_scale = brighter
    cutoff = min(freq * 8 * cutoff_scale, sr / 2 - 100)
    if cutoff > 100:
        sig = _lowpass(sig, cutoff, sr)

    env = _adsr(n, sr, attack=0.003, decay=0.15, sustain=0.35, release=0.15)
    return sig * env


def _synth_flute(freq: float, duration: float, sr: int,
                 cutoff_scale: float = 1.0) -> np.ndarray:
    """Flute: sine with vibrato LFO + filtered breath noise."""
    n = int(sr * duration)
    t = np.arange(n) / sr

    # Vibrato: 5Hz, increasing depth over time
    vibrato_depth = 0.003 * np.minimum(t / 0.3, 1.0)  # Ramp up over 300ms
    vibrato = vibrato_depth * np.sin(2 * np.pi * 5.0 * t)
    sig = np.sin(2 * np.pi * freq * (1.0 + vibrato) * t) * 0.8

    # Filtered breath noise with shaped envelope
    rng = np.random.default_rng(42)
    noise = rng.normal(0, 0.05, n)
    breath_env = _adsr(n, sr, attack=0.04, decay=0.1, sustain=0.3, release=0.1)
    cutoff = min(2000 * cutoff_scale, sr / 2 - 100)
    if cutoff > 100:
        noise = _lowpass(noise, cutoff, sr)
    sig += noise * breath_env

    env = _adsr(n, sr, attack=0.05, decay=0.05, sustain=0.6, release=0.12)
    return sig * env


def _synth_brass(freq: float, duration: float, sr: int,
                 cutoff_scale: float = 1.0) -> np.ndarray:
    """Brass: sawtooth through resonant LP with filter envelope."""
    n = int(sr * duration)
    t = np.arange(n) / sr

    sig = _sawtooth(freq, duration, sr) * 0.5

    # Filter envelope: opens on attack, then settles
    filter_env = _adsr(n, sr, attack=0.05, decay=0.2, sustain=0.4, release=0.1)
    # Sweep cutoff from low to high (stateful filtering for continuity)
    base_cutoff = freq * 2
    max_cutoff = min(freq * 10 * cutoff_scale, sr / 2 - 100)
    if max_cutoff > 100:
        block_size = sr // 20  # 50ms blocks
        filtered = np.zeros(n)
        zi = None
        for start in range(0, n, block_size):
            end = min(start + block_size, n)
            env_val = np.mean(filter_env[start:end])
            cutoff = base_cutoff + (max_cutoff - base_cutoff) * env_val
            cutoff = max(100, min(cutoff, sr / 2 - 100))
            filtered[start:end], zi = _lowpass_stateful(sig[start:end], cutoff, sr, zi=zi)
        sig = filtered

    env = _adsr(n, sr, attack=0.03, decay=0.1, sustain=0.5, release=0.08)
    return sig * env


def _synth_bass(freq: float, duration: float, sr: int,
                cutoff_scale: float = 1.0) -> np.ndarray:
    """Bass: filtered saw + sub sine + 2nd harmonic, with warm saturation."""
    n = int(sr * duration)
    t = np.arange(n) / sr

    # Band-limited sawtooth
    saw = _sawtooth(freq, duration, sr) * 0.30
    # Sub-oscillator at root frequency (not freq/2 which is often sub-audible)
    sub = _sine(freq, duration, sr) * 0.40
    # 2nd harmonic sine for body
    harm2 = _sine(freq * 2, duration, sr) * 0.15

    sig = saw + sub + harm2

    # Soft saturation for analog warmth
    _sat_drive = 1.5
    sig = np.tanh(sig * _sat_drive) / np.tanh(_sat_drive)

    # LP sweep: fast attack envelope (stateful filtering for continuity)
    filter_env = np.exp(-8.0 * t)
    base_cutoff = freq * 1.5
    max_cutoff = min(freq * 6 * cutoff_scale, sr / 2 - 100)
    if max_cutoff > 100:
        block_size = sr // 40  # 25ms blocks
        filtered = np.zeros(n)
        zi = None
        for start in range(0, n, block_size):
            end = min(start + block_size, n)
            env_val = np.mean(filter_env[start:end])
            cutoff = base_cutoff + (max_cutoff - base_cutoff) * env_val
            cutoff = max(100, min(cutoff, sr / 2 - 100))
            filtered[start:end], zi = _lowpass_stateful(
                sig[start:end], cutoff, sr, zi=zi, order=3)
        sig = filtered

    env = _adsr(n, sr, attack=0.005, decay=0.06, sustain=0.55, release=0.05)
    out = sig * env
    return out[:n]


def _synth_bell(freq: float, duration: float, sr: int,
                cutoff_scale: float = 1.0) -> np.ndarray:
    """Bell: FM synthesis with proper inharmonicity ratios."""
    n = int(sr * duration)
    t = np.arange(n) / sr

    # FM synthesis: carrier + modulator
    mod_ratio = 3.5  # Inharmonic ratio
    mod_index = 5.0 * np.exp(-3.0 * t)  # Decaying modulation index
    modulator = mod_index * np.sin(2 * np.pi * freq * mod_ratio * t)
    carrier = np.sin(2 * np.pi * freq * t + modulator)

    # Add secondary partials for richness
    sig = carrier * 0.6
    for ratio, amp, decay in [(2.76, 0.3, 2.0), (5.4, 0.15, 4.0), (8.93, 0.08, 6.0)]:
        sig += amp * np.sin(2 * np.pi * freq * ratio * t) * np.exp(-decay * t)

    env = _adsr(n, sr, attack=0.001, decay=0.4, sustain=0.05, release=0.4)
    return sig * env


def _synth_pad(freq: float, duration: float, sr: int,
               cutoff_scale: float = 1.0) -> np.ndarray:
    """Pad: 7 detuned saws (supersaw) through LP filter with slow LFO sweep."""
    n = int(sr * duration)
    t = np.arange(n) / sr

    detune_cents = [-12, -7, -3, 0, 3, 7, 12]
    sig = np.zeros(n)
    for dc in detune_cents:
        f = freq * (2.0 ** (dc / 1200.0))
        sig += _sawtooth(f, duration, sr) * (0.2 / len(detune_cents) * 2)

    # Slow LFO filter sweep (stateful filtering for continuity)
    lfo = 0.5 + 0.5 * np.sin(2 * np.pi * 0.2 * t)  # 0.2 Hz LFO
    base_cutoff = freq * 2
    max_cutoff = min(freq * 6 * cutoff_scale, sr / 2 - 100)
    if max_cutoff > 100:
        block_size = sr // 10  # 100ms blocks
        filtered = np.zeros(n)
        zi = None
        for start in range(0, n, block_size):
            end = min(start + block_size, n)
            lfo_val = np.mean(lfo[start:end])
            cutoff = base_cutoff + (max_cutoff - base_cutoff) * lfo_val
            cutoff = max(100, min(cutoff, sr / 2 - 100))
            filtered[start:end], zi = _lowpass_stateful(sig[start:end], cutoff, sr, zi=zi)
        sig = filtered

    env = _adsr(n, sr, attack=0.3, decay=0.15, sustain=0.75, release=0.4)
    return sig * env


def _synth_celesta(freq: float, duration: float, sr: int,
                   cutoff_scale: float = 1.0) -> np.ndarray:
    """Celesta: FM synthesis with fast decay, 2 modulators."""
    n = int(sr * duration)
    t = np.arange(n) / sr

    # Two FM modulators for bell-like timbre
    mod1_index = 3.0 * np.exp(-5.0 * t)
    mod2_index = 1.5 * np.exp(-8.0 * t)
    mod1 = mod1_index * np.sin(2 * np.pi * freq * 3.0 * t)
    mod2 = mod2_index * np.sin(2 * np.pi * freq * 5.0 * t)

    sig = np.sin(2 * np.pi * freq * t + mod1 + mod2)
    sig += 0.2 * np.sin(2 * np.pi * freq * 3 * t) * np.exp(-7.0 * t)

    env = _adsr(n, sr, attack=0.001, decay=0.2, sustain=0.15, release=0.2)
    return sig * env


def _synth_marimba(freq: float, duration: float, sr: int,
                   cutoff_scale: float = 1.0) -> np.ndarray:
    """Marimba: sine + subharmonic + short noise transient for mallet."""
    n = int(sr * duration)
    t = np.arange(n) / sr

    # Main tone with fast decay
    sig = np.sin(2 * np.pi * freq * t) * np.exp(-6.0 * t)
    # Subharmonic for body
    sig += 0.3 * np.sin(2 * np.pi * freq * 0.5 * t) * np.exp(-8.0 * t)
    # 4th harmonic
    sig += 0.4 * np.sin(2 * np.pi * freq * 4 * t) * np.exp(-15.0 * t)

    # Mallet transient (short filtered noise)
    rng = np.random.default_rng(42)
    click = rng.normal(0, 0.3, n)
    click_env = np.exp(-100.0 * t)  # Very fast decay
    cutoff = min(3000 * cutoff_scale, sr / 2 - 100)
    if cutoff > 100:
        click = _lowpass(click, cutoff, sr)
    sig += click * click_env

    return sig


def _synth_percussion(duration: float, sr: int,
                      cutoff_scale: float = 1.0) -> np.ndarray:
    """Percussion: layered sine body (pitch drop) + noise transient + filtered tail."""
    n = int(sr * duration)
    t = np.arange(n) / sr
    rng = np.random.default_rng(42)

    # Sine body with pitch drop (from 200Hz to 60Hz)
    freq_sweep = 200.0 * np.exp(-30.0 * t) + 60.0
    phase = np.cumsum(freq_sweep / sr) * 2 * np.pi
    body = np.sin(phase) * np.exp(-8.0 * t) * 0.5

    # Noise transient (short, bright)
    noise = rng.normal(0, 0.5, n)
    transient_env = np.exp(-50.0 * t)
    transient = noise * transient_env * 0.4

    # Filtered noise tail
    tail_noise = rng.normal(0, 0.3, n)
    cutoff = min(2000 * cutoff_scale, sr / 2 - 100)
    if cutoff > 100:
        tail_noise = _lowpass(tail_noise, cutoff, sr)
    tail_env = np.exp(-np.arange(n) / (sr * 0.08))
    tail = tail_noise * tail_env * 0.3

    return body + transient + tail


# ---------------------------------------------------------------------------
# Instrument routing
# ---------------------------------------------------------------------------

_INSTRUMENT_MAP: dict[str, callable] = {
    "acoustic_piano": _synth_piano,
    "synth_bass": _synth_bass,
    "counter_melody": _synth_flute,
    "synth_pad": _synth_pad,
    "motifs": _synth_bell,
}

_PROGRAM_MAP: dict[int, callable] = {
    0: _synth_piano,
    8: _synth_celesta,
    9: _synth_bell,       # Glockenspiel
    11: _synth_bell,      # Vibraphone
    12: _synth_marimba,
    38: _synth_bass,
    60: _synth_brass,
    73: _synth_flute,
    89: _synth_pad,
}


def _get_synth(layer: Layer, audio_params: dict | None = None) -> callable:
    """Get the synthesizer function for a layer."""
    if layer.instrument in _INSTRUMENT_MAP:
        return _INSTRUMENT_MAP[layer.instrument]
    if layer.program in _PROGRAM_MAP:
        return _PROGRAM_MAP[layer.program]
    return _synth_piano


# ---------------------------------------------------------------------------
# Layer synthesis (stereo + humanization)
# ---------------------------------------------------------------------------

def synthesize_layer(layer: Layer, feat: MolecularFeatures,
                     cfg: Config) -> np.ndarray:
    """Synthesize a single layer to a stereo numpy array (N, 2).

    Includes per-note panning, humanization (timing jitter, velocity variation),
    vibrato/pitch-bend effects, and filter warmth from audio parameters.
    """
    sr = cfg.audio_sample_rate
    ap = feat.audio_parameters
    bpm = ap["bpm"]
    duration_sec = ap["duration"]
    total_samples = int(sr * duration_sec) + sr  # Extra second for tails

    left = np.zeros(total_samples, dtype=np.float64)
    right = np.zeros(total_samples, dtype=np.float64)

    synth_fn = _get_synth(layer, ap)
    is_perc = layer.channel == 9

    # Humanization parameters
    swing = ap.get("swing", 0.0)
    filter_warmth = ap.get("filter_warmth", 0.5)
    cutoff_scale = 0.3 + 0.7 * (1.0 - filter_warmth)  # Warm = lower cutoff

    # Seeded RNG for reproducible humanization
    seed = ap.get("seed")
    rng = np.random.default_rng(seed if seed is not None else 42)

    for note in layer.notes:
        beat_time = note.start
        time_sec = beat_time * 60.0 / bpm
        start_sample = int(time_sec * sr)
        dur_sec = note.duration * 60.0 / bpm

        if start_sample >= total_samples or dur_sec <= 0:
            continue

        # Humanization: timing jitter (±10ms scaled by swing)
        jitter_samples = int(rng.uniform(-0.010, 0.010) * sr * (swing / 0.3 if swing > 0 else 0.1))
        start_sample = max(0, start_sample + jitter_samples)

        # Humanization: velocity variation (±5%)
        vel_variation = rng.uniform(-0.05, 0.05)
        velocity = max(1, min(127, int(note.velocity * (1.0 + vel_variation))))

        # Synthesize
        if is_perc:
            sig = _synth_percussion(dur_sec, sr, cutoff_scale=cutoff_scale)
        else:
            freq = _midi_to_freq(note.pitch)
            sig = synth_fn(freq, dur_sec, sr, cutoff_scale=cutoff_scale)

        # Per-note effects
        effects = note.effects
        if effects.get("vibrato"):
            from chemsymphony.synthesis.effects import apply_vibrato
            sig = apply_vibrato(sig, sr, rate=5.0, depth=0.003)

        if effects.get("pitch_bend"):
            # Slight pitch bend up: stretch signal slightly
            bend_amount = 1.02  # 2% up
            n_bent = int(len(sig) / bend_amount)
            if n_bent > 0:
                fractional_indices = np.linspace(0, len(sig) - 1, n_bent)
                sig = np.interp(fractional_indices, np.arange(len(sig)), sig)

        # Micro-fades (2ms cosine) to prevent clicks at note boundaries
        fade_samples = min(int(0.002 * sr), len(sig) // 2)
        if fade_samples > 0:
            fade_in = 0.5 * (1.0 - np.cos(np.linspace(0, np.pi, fade_samples)))
            fade_out = 0.5 * (1.0 + np.cos(np.linspace(0, np.pi, fade_samples)))
            sig[:fade_samples] *= fade_in
            sig[-fade_samples:] *= fade_out

        # Scale by velocity
        sig = sig * (velocity / 127.0)

        # Apply constant-power panning
        pan = max(-1.0, min(1.0, note.pan))
        angle = (pan + 1.0) * np.pi / 4.0
        l_gain = np.cos(angle)
        r_gain = np.sin(angle)

        # Place into stereo buffer
        end_sample = min(start_sample + len(sig), total_samples)
        length = end_sample - start_sample
        if length > 0 and start_sample >= 0:
            left[start_sample:end_sample] += sig[:length] * l_gain
            right[start_sample:end_sample] += sig[:length] * r_gain

    return np.column_stack([left, right])


def synthesize_audio(
    layers: CompositionLayers,
    feat: MolecularFeatures,
    cfg: Config,
) -> list[tuple[str, np.ndarray, float]]:
    """Synthesize all layers to stereo numpy arrays.

    Returns list of (layer_name, stereo_audio_array, mix_volume).
    """
    mix_cfg = cfg.mix
    ap = feat.audio_parameters
    arrangement_density = ap.get("arrangement_density", 1.0)

    result = []

    volume_map = {
        "lead_melody": mix_cfg.get("lead_melody", 0.85),
        "bass_loop": mix_cfg.get("bass_loops", 0.70),
        "counter_melody": mix_cfg.get("counter_melodies", 0.55),
        "drone_pad": mix_cfg.get("drone_pads", 0.40),
        "percussion": mix_cfg.get("percussion", 0.50),
        "motifs": mix_cfg.get("motifs", 0.65),
    }

    # Molecular-adaptive emphasis factors (0.7x–1.0x)
    # Emphasize layers that correspond to the molecule's dominant features
    chain_emphasis = min(1.0, 0.7 + 0.3 * (feat.longest_chain / max(feat.heavy_atom_count, 1)))
    ring_emphasis = min(1.0, 0.7 + 0.3 * (feat.ring_count / max(feat.heavy_atom_count / 3, 1)))
    aromatic_emphasis = min(1.0, 0.7 + 0.3 * feat.aromatic_fraction)
    branch_emphasis = min(1.0, 0.7 + 0.3 * (feat.branch_count / max(feat.heavy_atom_count / 3, 1)))
    fg_emphasis = min(1.0, 0.7 + 0.3 * (len(feat.functional_groups) / 5.0))

    emphasis_map = {
        "lead_melody": chain_emphasis,
        "bass_loop": ring_emphasis,
        "counter_melody": branch_emphasis,
        "drone_pad": aromatic_emphasis,
        "percussion": 1.0,  # Percussion always at full
        "motifs": fg_emphasis,
    }

    for layer in layers.all_layers():
        if not layer.notes:
            continue
        audio = synthesize_layer(layer, feat, cfg)

        # Determine mix volume from layer name
        vol = 0.5
        emphasis = 1.0
        for prefix, v in volume_map.items():
            if layer.name.startswith(prefix):
                vol = v
                emphasis = emphasis_map.get(prefix, 1.0)
                break

        # Scale by arrangement density and molecular emphasis
        vol *= arrangement_density * emphasis

        result.append((layer.name, audio, vol))

    return result
