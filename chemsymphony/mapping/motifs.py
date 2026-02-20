"""§8 mapping: Functional groups → musical motifs & accents."""

from __future__ import annotations

from chemsymphony.config import Config
from chemsymphony.features import MolecularFeatures
from chemsymphony.mapping.melody import Layer, NoteEvent

# Pre-designed motifs: each returns a list of (pitch_offset, start_offset, duration, velocity)
_MOTIF_LIBRARY: dict[str, list[tuple[int, float, float, int]]] = {
    # ── Alcohols & ethers ────────────────────────────────────────────
    "hydroxyl": [
        (0, 0.0, 0.25, 90),    # Root
        (4, 0.25, 0.25, 95),   # Major third up — bright arpeggio
    ],
    "phenol": [
        (0, 0.0, 0.2, 95),    # Bright arpeggio with resonant ring
        (4, 0.2, 0.2, 100),
        (7, 0.4, 0.6, 90),    # Sustained fifth — aromatic resonance
    ],
    "enol": [
        (0, 0.0, 0.4, 75),    # Wavering tone — tautomeric instability
        (1, 0.4, 0.4, 70),    # Semitone bend
        (0, 0.8, 0.4, 65),
    ],
    "ether": [
        (0, 0.0, 1.5, 60),    # Breathy sustained note
    ],
    "epoxide": [
        (0, 0.0, 0.1, 110),   # Tight burst — ring strain
        (5, 0.1, 0.1, 105),   # Tritone tension
        (0, 0.2, 0.1, 100),
    ],
    "peroxide": [
        (0, 0.0, 0.1, 115),   # Crackling double hit — reactive
        (0, 0.15, 0.1, 110),
        (12, 0.3, 0.1, 105),  # Octave spark
    ],
    "acetal": [
        (0, 0.0, 0.8, 60),    # Muffled bell — protected carbonyl
        (7, 0.0, 0.8, 50),    # Soft fifth underneath
    ],

    # ── Carbonyls ────────────────────────────────────────────────────
    "carbonyl": [
        (0, 0.0, 1.5, 80),    # Bell strike — long decay
    ],
    "aldehyde": [
        (0, 0.0, 1.2, 85),    # Bell + grace note
        (2, -0.1, 0.1, 70),
    ],
    "ketone": [
        (0, 0.0, 1.0, 75),    # Muted bell
    ],

    # ── Carboxylic acid & derivatives ────────────────────────────────
    "carboxylic_acid": [
        (0, 0.0, 1.0, 85),    # Bell strike
        (4, 1.0, 0.25, 90),   # + arpeggio
        (7, 1.25, 0.25, 90),
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
    "acyl_halide": [
        (0, 0.0, 0.15, 120),  # Sharp metallic strike — highly reactive
        (12, 0.15, 0.2, 110),  # Octave ping
    ],
    "anhydride": [
        (0, 0.0, 1.0, 80),    # Double bell — two carbonyls
        (0, 0.5, 1.0, 75),    # Second bell offset
        (7, 0.0, 1.5, 55),
    ],
    "carbamate": [
        (7, 0.0, 0.5, 70),    # Ester-like descending sigh
        (0, 0.5, 1.0, 65),    # + amine warmth sustained
        (4, 0.5, 1.0, 55),
    ],
    "urea": [
        (0, 0.0, 1.2, 65),    # Gentle double chord — two NH groups
        (4, 0.0, 1.2, 60),
        (3, 0.6, 1.2, 60),    # Minor second chord enters late
        (7, 0.6, 1.2, 55),
    ],
    "lactone": [
        (7, 0.0, 0.3, 80),    # Looping descending — cyclic ester
        (4, 0.3, 0.3, 75),
        (0, 0.6, 0.3, 70),
        (7, 0.9, 0.3, 75),    # Loop repeat
    ],
    "lactam": [
        (0, 0.0, 0.5, 75),    # Cyclic bell + pad loop
        (4, 0.0, 0.5, 60),
        (0, 0.5, 0.5, 70),    # Loop repeat
        (3, 0.5, 0.5, 55),
    ],

    # ── Nitrogen-containing ──────────────────────────────────────────
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
    "nitrile": [
        (24, 0.0, 2.0, 65),   # High sustained whine — 2 octaves up
    ],
    "isocyanate": [
        (12, 0.0, 0.3, 80),   # Rising whine — reactive N=C=O
        (18, 0.3, 0.3, 85),
        (24, 0.6, 1.0, 75),   # Sustained high
    ],
    "isothiocyanate": [
        (5, 0.0, 0.3, 75),    # Deeper rising whine — heavier S atom
        (12, 0.3, 0.3, 80),
        (17, 0.6, 1.0, 70),
    ],
    "nitro": [
        (0, 0.0, 0.15, 110),  # Sharp staccato burst
        (0, 0.2, 0.15, 105),
    ],
    "imine": [
        (0, 0.0, 0.5, 75),    # Upward slide — C=N bond
        (5, 0.3, 0.6, 80),    # Perfect fourth up
    ],
    "oxime": [
        (0, 0.0, 0.4, 75),    # Slide + bright ping
        (5, 0.3, 0.4, 80),
        (12, 0.7, 0.25, 90),  # OH sparkle
    ],
    "hydrazine": [
        (0, 0.0, 0.2, 85),    # Stuttering double — N-N
        (0, 0.25, 0.2, 80),
        (2, 0.5, 0.2, 80),    # Semitone drift
        (2, 0.75, 0.2, 75),
    ],
    "hydrazone": [
        (0, 0.0, 0.5, 75),    # Imine slide + stutter
        (5, 0.3, 0.4, 80),
        (5, 0.8, 0.2, 75),    # N-N stutter
        (5, 1.05, 0.2, 70),
    ],
    "azide": [
        (0, 0.0, 0.1, 115),   # Rapid triple hit — explosive N3
        (7, 0.12, 0.1, 110),
        (14, 0.24, 0.1, 105),
    ],
    "azo": [
        (0, 0.0, 1.5, 70),    # Buzzing sustained — N=N dye-like
        (1, 0.0, 1.5, 65),    # Semitone cluster = buzz
    ],
    "guanidine": [
        (0, 0.0, 1.2, 75),    # Rich major chord — strongly basic
        (4, 0.0, 1.2, 70),
        (7, 0.0, 1.2, 70),
        (12, 0.0, 1.2, 65),   # Octave brightness
    ],
    "enamine": [
        (0, 0.0, 0.6, 70),    # Smooth glide — N on C=C
        (2, 0.3, 0.6, 70),
        (5, 0.6, 0.6, 65),
    ],

    # ── Sulfur-containing ────────────────────────────────────────────
    "thiol": [
        (0, 0.0, 0.25, 110),  # Low growling accent
    ],
    "thioether": [
        (-5, 0.0, 1.2, 65),   # Dark sustained — S between carbons
    ],
    "disulfide": [
        (0, 0.0, 0.25, 105),  # Two low growls — S-S bridge
        (-5, 0.3, 0.25, 100),
    ],
    "sulfoxide": [
        (0, 0.0, 1.5, 70),    # Distorted pad swell
    ],
    "sulfone": [
        (0, 0.0, 1.5, 80),    # Heavy pad — fully oxidised S
        (-7, 0.0, 1.5, 70),   # Low fifth for weight
    ],
    "sulfonic_acid": [
        (0, 0.0, 0.2, 110),   # Bright aggressive hit — strong acid
        (4, 0.2, 0.2, 105),
        (7, 0.4, 0.8, 95),
    ],
    "sulfonamide": [
        (0, 0.0, 1.2, 75),    # Pad + warm chord — S(=O)2-N
        (3, 0.0, 1.2, 65),
        (7, 0.0, 1.2, 60),
    ],
    "sulfonate_ester": [
        (0, 0.0, 1.2, 75),    # Pad + descending — S(=O)2-O-C
        (7, 0.0, 0.5, 70),
        (0, 0.5, 0.5, 65),
    ],
    "thioamide": [
        (0, 0.0, 1.0, 70),    # Dark warm chord — C(=S)-N
        (3, 0.0, 1.0, 65),
        (-5, 0.0, 1.0, 60),
    ],
    "thioketone": [
        (-2, 0.0, 1.0, 70),   # Dark muted bell — C(=S), deeper than ketone
    ],

    # ── Phosphorus-containing ────────────────────────────────────────
    "phosphate": [
        (0, 0.0, 0.1, 85),    # Marimba roll
        (0, 0.1, 0.1, 80),
        (0, 0.2, 0.1, 75),
        (0, 0.3, 0.1, 70),
    ],
    "phosphonate": [
        (0, 0.0, 0.1, 80),    # Marimba variant — with C attachment
        (5, 0.1, 0.1, 75),
        (0, 0.2, 0.1, 70),
        (5, 0.3, 0.1, 65),
    ],
    "phosphine": [
        (-7, 0.0, 1.5, 60),   # Low warm tone — P with 3 carbons
    ],
    "phosphoramide": [
        (0, 0.0, 0.1, 80),    # Marimba + chord — P with N
        (0, 0.1, 0.1, 75),
        (3, 0.2, 1.0, 65),
        (7, 0.2, 1.0, 60),
    ],

    # ── Boron-containing ─────────────────────────────────────────────
    "boronic_acid": [
        (12, 0.0, 0.3, 85),   # Glassy ping — B(OH)2
        (16, 0.3, 0.3, 80),
    ],
    "boronate_ester": [
        (16, 0.0, 0.4, 80),   # Glassy descending — B(OR)2
        (12, 0.4, 0.4, 75),
    ],

    # ── Halides ──────────────────────────────────────────────────────
    "halide": [
        (12, 0.0, 0.25, 95),  # Metallic ping — octave up
    ],

    # ── Unsaturated ──────────────────────────────────────────────────
    "alkene": [
        (0, 0.0, 0.4, 70),    # Smooth ascending slide — C=C
        (2, 0.2, 0.4, 70),
    ],
    "alkyne": [
        (24, 0.0, 1.5, 70),   # High piercing sustained — C≡C
        (19, 0.0, 1.5, 60),   # + fifth below for shimmer
    ],

    # ── Heterocyclic ─────────────────────────────────────────────────
    "pyridine": [
        (0, 0.0, 0.4, 70),    # Mellow 6-note loop — aromatic N ring
        (2, 0.4, 0.4, 65),
        (4, 0.8, 0.4, 65),
        (5, 1.2, 0.4, 60),
        (4, 1.6, 0.4, 60),
        (2, 2.0, 0.4, 55),
    ],
    "pyrrole": [
        (0, 0.0, 1.5, 60),    # Warm drone — NH in 5-ring
        (7, 0.0, 1.5, 50),    # Soft fifth
    ],
    "furan": [
        (12, 0.0, 0.2, 75),   # Airy flutter — O in 5-ring
        (11, 0.2, 0.2, 70),
        (12, 0.4, 0.2, 65),
        (14, 0.6, 0.2, 70),
        (12, 0.8, 0.2, 60),
    ],
    "thiophene": [
        (-5, 0.0, 1.5, 65),   # Dark drone — S in 5-ring
        (-12, 0.0, 1.5, 55),  # Low octave
    ],
    "imidazole": [
        (0, 0.0, 0.8, 70),    # Complex chord — two N in ring
        (3, 0.0, 0.8, 65),
        (7, 0.0, 0.8, 65),
        (10, 0.0, 0.8, 60),
    ],
    "oxazole": [
        (12, 0.0, 0.3, 75),   # Airy + bright — O and N
        (7, 0.3, 0.5, 70),
        (4, 0.8, 0.3, 70),
    ],
    "thiazole": [
        (-5, 0.0, 0.5, 70),   # Dark + bright — S and N
        (7, 0.5, 0.5, 75),
        (0, 1.0, 0.5, 65),
    ],
    "indole": [
        (0, 0.0, 2.0, 60),    # Rich warm pad — fused pyrrole+benzene
        (4, 0.0, 2.0, 55),
        (7, 0.0, 2.0, 55),
        (12, 0.0, 2.0, 50),
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
