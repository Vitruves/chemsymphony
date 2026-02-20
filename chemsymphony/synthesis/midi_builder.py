"""MIDI track assembly using midiutil."""

from __future__ import annotations

import io
from typing import Any

from midiutil import MIDIFile

from chemsymphony.config import Config
from chemsymphony.features import MolecularFeatures
from chemsymphony.mapping import CompositionLayers
from chemsymphony.mapping.melody import Layer


def _layer_to_midi(layer: Layer, bpm: int) -> bytes:
    """Convert a single Layer to MIDI bytes."""
    midi = MIDIFile(1, deinterleave=False)
    midi.addTempo(0, 0, bpm)

    # Set instrument
    channel = layer.channel
    if channel != 9:  # Don't set program on percussion channel
        midi.addProgramChange(0, channel, 0, layer.program)

    for note in layer.notes:
        midi.addNote(
            track=0,
            channel=channel,
            pitch=note.pitch,
            time=note.start,
            duration=max(0.01, note.duration),
            volume=note.velocity,
        )

    buf = io.BytesIO()
    midi.writeFile(buf)
    return buf.getvalue()


def build_midi(
    layers: CompositionLayers,
    feat: MolecularFeatures,
    cfg: Config,
) -> dict[str, Any]:
    """Build MIDI data for all layers.

    Returns a dict with:
      - "master": bytes (combined MIDI with all tracks)
      - "tracks": list of (name, bytes) tuples (individual MIDI files)
    """
    ap = feat.audio_parameters
    bpm = ap["bpm"]
    all_layers = layers.all_layers()

    # Build individual track MIDI files
    track_list: list[tuple[str, bytes]] = []
    for layer in all_layers:
        if layer.notes:
            data = _layer_to_midi(layer, bpm)
            track_list.append((layer.name, data))

    # Build master MIDI with all tracks
    n_tracks = len(all_layers)
    if n_tracks == 0:
        n_tracks = 1
    master = MIDIFile(n_tracks, deinterleave=False)
    master.addTempo(0, 0, bpm)

    for track_idx, layer in enumerate(all_layers):
        channel = layer.channel
        if channel != 9:
            master.addProgramChange(track_idx, channel, 0, layer.program)

        for note in layer.notes:
            master.addNote(
                track=track_idx,
                channel=channel,
                pitch=note.pitch,
                time=note.start,
                duration=max(0.01, note.duration),
                volume=note.velocity,
            )

    buf = io.BytesIO()
    master.writeFile(buf)
    master_bytes = buf.getvalue()

    return {
        "master": master_bytes,
        "tracks": track_list,
    }
