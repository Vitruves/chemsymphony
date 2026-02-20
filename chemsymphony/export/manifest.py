"""JSON manifest generation."""

from __future__ import annotations

from typing import Any

from chemsymphony.canonicalize import CanonicalResult
from chemsymphony.config import Config
from chemsymphony.features import MolecularFeatures
from chemsymphony.mapping import CompositionLayers


def build_manifest(
    canon: CanonicalResult,
    feat: MolecularFeatures,
    layers: CompositionLayers,
    cfg: Config,
) -> dict[str, Any]:
    """Build the manifest dict for MIDI export."""
    ap = feat.audio_parameters

    tracks_info = []
    for i, layer in enumerate(layers.all_layers(), start=1):
        if not layer.notes:
            continue
        info: dict[str, Any] = {
            "name": layer.name,
            "file": f"tracks/{i:02d}_{layer.name}.mid",
            "channel": layer.channel,
            "instrument": layer.instrument,
            "notes": len(layer.notes),
        }
        # Add bass loop specific info
        if layer.name.startswith("bass_loop"):
            ring_idx = int(layer.name.split("ring")[-1]) - 1 if "ring" in layer.name else 0
            if ring_idx < len(feat.rings):
                ring = feat.rings[ring_idx]
                info["loop_beats"] = ring["size"]
                info["aromatic"] = ring["is_aromatic"]
                info["source_feature"] = f"ring_system_{ring_idx}"
        elif layer.name == "lead_melody":
            info["source_feature"] = "main_carbon_chain"

        tracks_info.append(info)

    return {
        "smiles_input": canon.input_smiles,
        "smiles_canonical": canon.canonical_smiles,
        "molecular_formula": feat.molecular_formula,
        "features_extracted": {
            "heavy_atom_count": feat.heavy_atom_count,
            "molecular_weight": round(feat.molecular_weight, 2),
            "ring_count": feat.ring_count,
            "aromatic_ring_count": feat.aromatic_ring_count,
            "longest_chain": feat.longest_chain,
            "branch_count": feat.branch_count,
            "stereo_centers": feat.stereo_center_count,
            "functional_groups": list({fg["name"] for fg in feat.functional_groups}),
            "graph_diameter": feat.graph_diameter,
            "unique_elements": feat.unique_elements,
        },
        "audio_parameters": {
            "bpm": ap["bpm"],
            "key": ap["root_note"],
            "scale": ap["scale"],
            "duration_seconds": round(ap["duration"], 2),
            "root_note_midi": ap["root_note_midi"],
            "harmonic_density": ap["harmonic_density"],
            "note_density": round(ap["note_density"], 2),
        },
        "tracks": tracks_info,
    }
