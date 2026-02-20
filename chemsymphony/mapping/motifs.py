"""§8 mapping: Functional groups → musical motifs & accents."""

from __future__ import annotations

from chemsymphony.config import Config
from chemsymphony.features import MolecularFeatures
from chemsymphony.mapping.melody import Layer, NoteEvent

# Pre-designed motifs: each returns a list of (pitch_offset, start_offset, duration, velocity)
_MOTIF_LIBRARY: dict[str, list[tuple[int, float, float, int]]] = {
    "hydroxyl": [
        (0, 0.0, 0.25, 90),    # Root
        (4, 0.25, 0.25, 95),   # Major third up — bright arpeggio
    ],
    "carbonyl": [
        (0, 0.0, 1.5, 80),    # Bell strike — long decay
    ],
    "carboxylic_acid": [
        (0, 0.0, 1.0, 85),    # Bell strike
        (4, 1.0, 0.25, 90),   # + arpeggio
        (7, 1.25, 0.25, 90),
    ],
    "primary_amine": [
        (0, 0.0, 1.0, 70),    # Warm chord swell — major
        (4, 0.0, 1.0, 65),
        (7, 0.0, 1.0, 60),
    ],
    "secondary_amine": [
        (0, 0.0, 1.0, 70),    # Minor chord
        (3, 0.0, 1.0, 65),
        (7, 0.0, 1.0, 60),
    ],
    "tertiary_amine": [
        (0, 0.0, 1.0, 70),    # Diminished chord
        (3, 0.0, 1.0, 65),
        (6, 0.0, 1.0, 60),
    ],
    "ester": [
        (7, 0.0, 0.5, 75),    # Descending two-note figure
        (0, 0.5, 0.5, 70),
    ],
    "amide": [
        (0, 0.0, 1.0, 75),    # Bell + pad swell
        (4, 0.0, 1.5, 60),
        (7, 0.0, 1.5, 55),
    ],
    "ether": [
        (0, 0.0, 1.5, 60),    # Breathy sustained note
    ],
    "thiol": [
        (0, 0.0, 0.25, 110),  # Low growling accent
    ],
    "nitro": [
        (0, 0.0, 0.15, 110),  # Sharp staccato burst
        (0, 0.2, 0.15, 105),
    ],
    "phosphate": [
        (0, 0.0, 0.1, 85),    # Marimba roll
        (0, 0.1, 0.1, 80),
        (0, 0.2, 0.1, 75),
        (0, 0.3, 0.1, 70),
    ],
    "halide": [
        (12, 0.0, 0.25, 95),  # Metallic ping — octave up
    ],
    "aldehyde": [
        (0, 0.0, 1.2, 85),    # Bell + grace note
        (2, -0.1, 0.1, 70),
    ],
    "ketone": [
        (0, 0.0, 1.0, 75),    # Muted bell
    ],
    "nitrile": [
        (24, 0.0, 2.0, 65),   # High sustained whine — 2 octaves up
    ],
    "sulfoxide": [
        (0, 0.0, 1.5, 70),    # Distorted pad swell
    ],
}


def generate_motifs(feat: MolecularFeatures, cfg: Config) -> Layer:
    """Generate functional group motif layer."""
    ap = feat.audio_parameters
    root = ap["root_note_midi"]
    bpm = ap["bpm"]
    duration_sec = ap["duration"]
    beats_total = (duration_sec / 60.0) * bpm
    chain_len = max(feat.longest_chain, 1)

    notes: list[NoteEvent] = []

    for fg in feat.functional_groups:
        name = fg["name"]
        atoms = fg["atoms"]
        motif = _MOTIF_LIBRARY.get(name)
        if motif is None:
            continue

        # Position: earliest atom in the functional group, mapped to time
        if atoms and feat.longest_chain_atoms:
            # Find the chain position of the closest atom
            chain_set = set(feat.longest_chain_atoms)
            min_pos = 0
            for a in atoms:
                if a in chain_set:
                    min_pos = feat.longest_chain_atoms.index(a)
                    break
            entry_beat = (min_pos / chain_len) * beats_total
        else:
            entry_beat = 0.0

        for pitch_off, start_off, dur, vel in motif:
            t = max(0.0, entry_beat + start_off)
            notes.append(NoteEvent(
                pitch=max(36, min(108, root + pitch_off)),
                start=t,
                duration=dur,
                velocity=vel,
                channel=5,
            ))

    return Layer(
        name="motifs",
        instrument="motifs",
        channel=5,
        notes=notes,
        program=11,  # Vibraphone (close to bell)
    )
