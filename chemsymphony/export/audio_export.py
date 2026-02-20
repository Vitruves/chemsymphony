"""MP3/WAV export."""

from __future__ import annotations

import io
import struct

import numpy as np

from chemsymphony.config import Config


def _to_wav_bytes(audio: np.ndarray, sr: int, bit_depth: int = 16) -> bytes:
    """Convert a numpy stereo array to WAV file bytes."""
    n_channels = audio.shape[1] if audio.ndim == 2 else 1
    n_frames = audio.shape[0]

    # Scale to int range
    if bit_depth == 16:
        max_val = 32767
        dtype = np.int16
    else:
        max_val = 32767
        dtype = np.int16

    scaled = np.clip(audio * max_val, -max_val, max_val).astype(dtype)
    raw_data = scaled.tobytes()

    # Build WAV header
    bytes_per_sample = bit_depth // 8
    block_align = n_channels * bytes_per_sample
    byte_rate = sr * block_align
    data_size = len(raw_data)
    file_size = 36 + data_size

    buf = io.BytesIO()
    buf.write(b"RIFF")
    buf.write(struct.pack("<I", file_size))
    buf.write(b"WAVE")
    buf.write(b"fmt ")
    buf.write(struct.pack("<I", 16))          # Chunk size
    buf.write(struct.pack("<H", 1))           # PCM format
    buf.write(struct.pack("<H", n_channels))
    buf.write(struct.pack("<I", sr))
    buf.write(struct.pack("<I", byte_rate))
    buf.write(struct.pack("<H", block_align))
    buf.write(struct.pack("<H", bit_depth))
    buf.write(b"data")
    buf.write(struct.pack("<I", data_size))
    buf.write(raw_data)

    return buf.getvalue()


def export_audio(audio: np.ndarray, *, fmt: str = "mp3", config: Config) -> bytes:
    """Export a stereo numpy array to audio bytes (WAV or MP3)."""
    sr = config.audio_sample_rate
    bit_depth = config.audio_bit_depth

    wav_bytes = _to_wav_bytes(audio, sr, bit_depth)

    if fmt == "wav":
        return wav_bytes

    if fmt == "mp3":
        try:
            from pydub import AudioSegment
            segment = AudioSegment.from_wav(io.BytesIO(wav_bytes))
            mp3_buf = io.BytesIO()
            segment.export(mp3_buf, format="mp3", bitrate=f"{config.audio_mp3_bitrate}k")
            return mp3_buf.getvalue()
        except Exception:
            # Fallback: return WAV if pydub/ffmpeg unavailable
            return wav_bytes

    return wav_bytes
