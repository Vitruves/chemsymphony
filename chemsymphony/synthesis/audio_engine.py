"""Direct audio synthesis engine using numpy/scipy."""

from __future__ import annotations

import numpy as np

from chemsymphony.config import Config
from chemsymphony.features import MolecularFeatures
from chemsymphony.mapping import CompositionLayers
from chemsymphony.mapping.melody import Layer, NoteEvent

# --- Oscillators ---

def _sine(freq: float, duration: float, sr: int) -> np.ndarray:
    t = np.arange(int(sr * duration)) / sr
    return np.sin(2 * np.pi * freq * t)


def _square(freq: float, duration: float, sr: int) -> np.ndarray:
    return np.sign(_sine(freq, duration, sr))


def _sawtooth(freq: float, duration: float, sr: int) -> np.ndarray:
    t = np.arange(int(sr * duration)) / sr
    return 2.0 * (freq * t % 1.0) - 1.0


def _triangle(freq: float, duration: float, sr: int) -> np.ndarray:
    return 2.0 * np.abs(_sawtooth(freq, duration, sr)) - 1.0


# --- ADSR Envelope ---

def _adsr(n_samples: int, sr: int, attack: float = 0.01, decay: float = 0.05,
          sustain: float = 0.7, release: float = 0.05) -> np.ndarray:
    """Generate ADSR envelope, clamped to fit within n_samples."""
    if n_samples <= 0:
        return np.zeros(0)

    a = int(attack * sr)
    d = int(decay * sr)
    r = int(release * sr)

    # Clamp so total A+D+R doesn't exceed n_samples
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
        env[pos:pos + d] = np.linspace(1, sustain, d)
        pos += d
    if s > 0:
        env[pos:pos + s] = sustain
        pos += s
    if r > 0:
        env[pos:pos + r] = np.linspace(sustain, 0, r)
    return env


# --- Instruments ---

def _midi_to_freq(midi_note: int) -> float:
    return 440.0 * (2.0 ** ((midi_note - 69) / 12.0))


def _synth_piano(freq: float, duration: float, sr: int) -> np.ndarray:
    """Additive synthesis piano: fundamental + harmonics."""
    n = int(sr * duration)
    sig = np.zeros(n)
    for k, amp in enumerate([1.0, 0.5, 0.3, 0.2, 0.1], start=1):
        t = np.arange(n) / sr
        sig += amp * np.sin(2 * np.pi * freq * k * t) * np.exp(-k * 0.5 * t)
    env = _adsr(n, sr, attack=0.005, decay=0.1, sustain=0.4, release=0.1)
    return sig * env


def _synth_flute(freq: float, duration: float, sr: int) -> np.ndarray:
    """Sine + breath noise."""
    n = int(sr * duration)
    sig = _sine(freq, duration, sr) * 0.8
    noise = np.random.default_rng(42).normal(0, 0.03, n)
    env = _adsr(n, sr, attack=0.05, decay=0.05, sustain=0.6, release=0.1)
    return (sig + noise) * env


def _synth_brass(freq: float, duration: float, sr: int) -> np.ndarray:
    """Sawtooth with filter-like envelope."""
    n = int(sr * duration)
    sig = _sawtooth(freq, duration, sr) * 0.5
    env = _adsr(n, sr, attack=0.03, decay=0.1, sustain=0.5, release=0.08)
    return sig * env


def _synth_bass(freq: float, duration: float, sr: int) -> np.ndarray:
    """Square + sub-octave."""
    n = int(sr * duration)
    sig = _square(freq, duration, sr) * 0.4
    sub = _sine(freq / 2, duration, sr) * 0.4
    env = _adsr(n, sr, attack=0.005, decay=0.05, sustain=0.6, release=0.05)
    out = (sig + sub) * env
    return out[:n]


def _synth_bell(freq: float, duration: float, sr: int) -> np.ndarray:
    """Inharmonic partials for bell tone."""
    n = int(sr * duration)
    sig = np.zeros(n)
    t = np.arange(n) / sr
    for ratio, amp in [(1.0, 1.0), (2.76, 0.5), (5.4, 0.3), (8.93, 0.15)]:
        sig += amp * np.sin(2 * np.pi * freq * ratio * t) * np.exp(-ratio * t)
    env = _adsr(n, sr, attack=0.001, decay=0.3, sustain=0.1, release=0.3)
    return sig * env


def _synth_pad(freq: float, duration: float, sr: int) -> np.ndarray:
    """Detuned sines for pad."""
    n = int(sr * duration)
    sig = np.zeros(n)
    for detune in [-3, 0, 3, 5]:
        f = freq * (2.0 ** (detune / 1200.0))
        sig += _sine(f, duration, sr) * 0.3
    env = _adsr(n, sr, attack=0.2, decay=0.1, sustain=0.8, release=0.3)
    return sig * env


def _synth_celesta(freq: float, duration: float, sr: int) -> np.ndarray:
    """Bright bell-like celesta."""
    n = int(sr * duration)
    t = np.arange(n) / sr
    sig = np.sin(2 * np.pi * freq * t) * np.exp(-3.0 * t)
    sig += 0.3 * np.sin(2 * np.pi * freq * 3 * t) * np.exp(-5.0 * t)
    env = _adsr(n, sr, attack=0.001, decay=0.15, sustain=0.2, release=0.15)
    return sig * env


def _synth_marimba(freq: float, duration: float, sr: int) -> np.ndarray:
    """Percussive pitched tone."""
    n = int(sr * duration)
    t = np.arange(n) / sr
    sig = np.sin(2 * np.pi * freq * t) * np.exp(-8.0 * t)
    sig += 0.5 * np.sin(2 * np.pi * freq * 4 * t) * np.exp(-12.0 * t)
    return sig


def _synth_percussion(duration: float, sr: int) -> np.ndarray:
    """Filtered noise burst for percussion."""
    n = int(sr * duration)
    rng = np.random.default_rng(42)
    noise = rng.normal(0, 0.5, n)
    env = np.exp(-np.arange(n) / (sr * 0.05))
    return noise * env


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
    9: _synth_bell,  # Glockenspiel
    11: _synth_bell,  # Vibraphone
    12: _synth_marimba,
    38: _synth_bass,
    60: _synth_brass,
    73: _synth_flute,
    89: _synth_pad,
}


def _get_synth(layer: Layer) -> callable:
    """Get the synthesizer function for a layer."""
    if layer.instrument in _INSTRUMENT_MAP:
        return _INSTRUMENT_MAP[layer.instrument]
    if layer.program in _PROGRAM_MAP:
        return _PROGRAM_MAP[layer.program]
    return _synth_piano


def synthesize_layer(layer: Layer, feat: MolecularFeatures, cfg: Config) -> np.ndarray:
    """Synthesize a single layer to a numpy audio array."""
    sr = cfg.audio_sample_rate
    ap = feat.audio_parameters
    bpm = ap["bpm"]
    duration_sec = ap["duration"]
    total_samples = int(sr * duration_sec) + sr  # Extra second for tails

    output = np.zeros(total_samples, dtype=np.float64)
    synth_fn = _get_synth(layer)

    is_perc = layer.channel == 9

    for note in layer.notes:
        beat_time = note.start
        time_sec = beat_time * 60.0 / bpm
        start_sample = int(time_sec * sr)
        dur_sec = note.duration * 60.0 / bpm

        if start_sample >= total_samples or dur_sec <= 0:
            continue

        if is_perc:
            sig = _synth_percussion(dur_sec, sr)
        else:
            freq = _midi_to_freq(note.pitch)
            sig = synth_fn(freq, dur_sec, sr)

        # Scale by velocity
        sig = sig * (note.velocity / 127.0)

        # Place into output buffer
        end_sample = min(start_sample + len(sig), total_samples)
        length = end_sample - start_sample
        if length > 0 and start_sample >= 0:
            output[start_sample:end_sample] += sig[:length]

    return output


def synthesize_audio(
    layers: CompositionLayers,
    feat: MolecularFeatures,
    cfg: Config,
) -> list[tuple[str, np.ndarray, float]]:
    """Synthesize all layers to numpy arrays.

    Returns list of (layer_name, audio_array, mix_volume).
    """
    mix_cfg = cfg.mix
    result = []

    volume_map = {
        "lead_melody": mix_cfg.get("lead_melody", 0.85),
        "bass_loop": mix_cfg.get("bass_loops", 0.70),
        "counter_melody": mix_cfg.get("counter_melodies", 0.55),
        "drone_pad": mix_cfg.get("drone_pads", 0.40),
        "percussion": mix_cfg.get("percussion", 0.50),
        "motifs": mix_cfg.get("motifs", 0.65),
    }

    for layer in layers.all_layers():
        if not layer.notes:
            continue
        audio = synthesize_layer(layer, feat, cfg)

        # Determine mix volume from layer name
        vol = 0.5
        for prefix, v in volume_map.items():
            if layer.name.startswith(prefix):
                vol = v
                break

        result.append((layer.name, audio, vol))

    return result
