<p align="center">
  <img src="chemsymphony.png" alt="ChemSymphony Logo" width="300">
</p>

# 🎵 ChemSymphony

**Generate unique audio from molecular SMILES strings.**

ChemSymphony transforms the structural information encoded in SMILES (Simplified Molecular-Input Line-Entry System) strings into rich, layered audio compositions. Simple molecules produce minimal melodies; complex molecules produce dense, intricate soundscapes. Every structural feature of a molecule — its atoms, bonds, rings, chains, branches, stereochemistry, and aromaticity — is mapped to a distinct musical dimension, ensuring that each molecule has its own sonic fingerprint.

---

## Table of Contents

- [Installation](#installation)
- [Quick Start](#quick-start)
- [CLI Reference](#cli-reference)
- [SMILES Canonicalization](#smiles-canonicalization)
- [Architecture Overview](#architecture-overview)
- [SMILES-to-Audio Mapping Pipeline](#smiles-to-audio-mapping-pipeline)
  - [1. Global Molecular Properties → Master Parameters](#1-global-molecular-properties--master-parameters)
  - [2. Atomic Composition → Instrument Palette & Timbre](#2-atomic-composition--instrument-palette--timbre)
  - [3. Ring Systems → Bass Loops & Rhythmic Layers](#3-ring-systems--bass-loops--rhythmic-layers)
  - [4. Aromaticity → Drone Pads & Reverb](#4-aromaticity--drone-pads--reverb)
  - [5. Main Carbon Chain → Lead Melody](#5-main-carbon-chain--lead-melody)
  - [6. Branches & Ramifications → Counter-Melodies & Harmonies](#6-branches--ramifications--counter-melodies--harmonies)
  - [7. Bond Types → Articulation & Dynamics](#7-bond-types--articulation--dynamics)
  - [8. Functional Groups → Motifs & Accents](#8-functional-groups--motifs--accents)
  - [9. Stereochemistry → Spatial Audio & Phrasing](#9-stereochemistry--spatial-audio--phrasing)
  - [10. Atom Connectivity & Graph Topology → Structure & Form](#10-atom-connectivity--graph-topology--structure--form)
  - [11. Atom Recurrence & Elemental Distribution → Rhythmic Patterns](#11-atom-recurrence--elemental-distribution--rhythmic-patterns)
  - [12. Charges & Electronegativity → Expression & Effects](#12-charges--electronegativity--expression--effects)
- [Output Formats](#output-formats)
- [Examples](#examples)
- [Python API](#python-api)
- [Configuration](#configuration)
- [Dependencies](#dependencies)
- [License](#license)

---

## Installation

ChemSymphony is a standard Python package installable via pip (or uv):

```bash
# With pip
pip install .

# With uv
uv pip install .

# Development install
pip install -e ".[dev]"
```

### Requirements

- Python ≥ 3.10
- RDKit (cheminformatics backbone)
- midiutil (MIDI generation)
- NumPy / SciPy (signal processing & audio synthesis)
- pydub + ffmpeg (MP3/WAV export)

---

## Quick Start

```bash
# Generate an MP3 from caffeine
chemsymphony -s "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"

# Output a WAV file to a specific path
chemsymphony -s "CCO" -o ethanol.wav -f wav

# Export MIDI tracks + JSON manifest
chemsymphony -s "c1ccccc1" -f midi -o benzene_tracks/

# Verbose mode — show all extracted features
chemsymphony -s "CC(=O)OC1=CC=CC=C1C(=O)O" -v
```

---

## CLI Reference

```
chemsymphony [options] -s/--smiles <SMILES>

Required:
  -s, --smiles TEXT          Input SMILES string

Output options:
  -o, --output PATH          Output file or directory (default: ./output.<format>)
  -f, --format FORMAT        Output format: mp3 | wav | midi (default: mp3)

Audio tuning:
  --bpm INT                  Override base tempo (default: auto-derived from molecule)
  --duration FLOAT           Override duration in seconds (default: auto-derived)
  --key TEXT                 Force musical key, e.g. "Cm", "G" (default: auto-derived)
  --seed INT                 Random seed for reproducible stochastic elements

Verbosity:
  -v, --verbose              Print extracted molecular features and mapping details
  -q, --quiet                Suppress all non-error output
  --dry-run                  Extract features and print mapping without generating audio

General:
  --version                  Show version and exit
  -h, --help                 Show this help message and exit
```

---

## SMILES Canonicalization

**The very first step in the pipeline is SMILES canonicalization.** Before any structural analysis begins, the input SMILES string is converted to its canonical form using RDKit's `Chem.MolToSmiles(Chem.MolFromSmiles(input), canonical=True)`. This guarantees that any valid SMILES representation of the same molecule produces identical audio output:

```
"OCC"  →  canonical: "CCO"  →  audio X
"C(O)C" →  canonical: "CCO"  →  audio X   (identical)
```

Canonicalization also serves as input validation — invalid SMILES strings are rejected with a clear error message before any processing occurs.

---

## Architecture Overview

The pipeline flows through four stages:

```
┌──────────────┐    ┌───────────────────┐    ┌─────────────────┐    ┌──────────────┐
│  SMILES      │───▶│  Feature          │───▶│  Audio          │───▶│  Renderer    │
│  Canonicalize│    │  Extraction       │    │  Mapping Engine │    │  & Export    │
│  & Parse     │    │  (12 pipelines)   │    │  (per-layer     │    │  mp3/wav/midi│
└──────────────┘    └───────────────────┘    │   composition)  │    └──────────────┘
                                             └─────────────────┘
```

1. **Canonicalize & Parse** — Normalize SMILES, build the molecular graph via RDKit.
2. **Feature Extraction** — Twelve parallel analysis pipelines each extract a specific class of structural information (atoms, rings, chains, branches, bonds, etc.).
3. **Audio Mapping Engine** — Each feature class is mapped to a dedicated audio layer (melody, bass, pads, percussion, effects, etc.) with its own synthesis rules.
4. **Renderer & Export** — Layers are mixed and rendered to the requested output format.

---

## SMILES-to-Audio Mapping Pipeline

Below is an exhaustive description of every information-to-sound pipeline. Each section describes *what* structural information is extracted, *how* it maps to audio parameters, and *why* the mapping is musically meaningful.

---

### 1. Global Molecular Properties → Master Parameters

Global properties set the overarching musical parameters — the "world" in which all other layers operate.

| Molecular Feature | Audio Parameter | Mapping Logic |
|---|---|---|
| **Heavy atom count** | **Duration** | More atoms → longer composition. Linear scale: 2s baseline + 0.3s per heavy atom, capped at 120s. |
| **Molecular weight** | **Tempo (BPM)** | Light molecules (MW < 100): slow tempo (70–90 BPM). Medium (100–500): moderate (90–130 BPM). Heavy (>500): fast, dense (130–180 BPM). This creates simple/sparse feel for small molecules and urgency for large ones. |
| **Heteroatom-to-carbon ratio** | **Musical scale / mode** | Pure hydrocarbons (ratio ≈ 0): pentatonic major (bright, simple). Low heteroatom ratio (<0.3): major/Ionian. Moderate (0.3–0.6): Mixolydian or Dorian. High (>0.6): minor, Phrygian, or chromatic. Molecules with many heteroatoms sound darker and more dissonant. |
| **Atom type diversity** (unique element count) | **Harmonic density** | 1–2 element types: monophonic or simple harmony. 3–4: triads and simple chords. 5+: extended chords (7ths, 9ths), denser voicings. |
| **Total bond count** | **Note density / event rate** | More bonds → more MIDI events per beat. Sparse molecules breathe; dense ones are packed with notes. |
| **Molecular formula string hash** | **Root note** | The canonical molecular formula is hashed to select the root pitch (C, C#, D, … B), guaranteeing deterministic key assignment per molecule. |

**Musical rationale:** These master parameters create an immediate, intuitive relationship: small simple molecules (e.g., methane, ethanol) feel light and sparse; large complex molecules (e.g., taxol, chlorophyll) feel dense, fast, and harmonically rich.

---

### 2. Atomic Composition → Instrument Palette & Timbre

Each element present in the molecule contributes a dedicated instrument voice to the mix. The *quantity* of each element controls that voice's prominence (volume and note density).

| Element | Instrument / Timbre | Volume Scaling |
|---|---|---|
| **Carbon (C)** | Acoustic piano or clean electric piano — the "backbone" voice | Always present; volume proportional to C fraction |
| **Hydrogen (implicit)** | High-register soft chime or celesta — airy presence | Scaled by H count; very light in the mix |
| **Oxygen (O)** | Warm flute or woodwind — "breath" quality | Louder with more O atoms |
| **Nitrogen (N)** | Muted brass or French horn — resonant warmth | Scales with N count |
| **Sulfur (S)** | Distorted or growling synth bass — dark, heavy | Prominent when S is present |
| **Phosphorus (P)** | Marimba or xylophone — percussive and pitched | One voice per P atom |
| **Halogens (F, Cl, Br, I)** | Metallic percussion: hi-hat (F), snare (Cl), crash (Br), gong (I) | Each halogen adds a distinct percussive layer |
| **Metals (Na, K, Fe, etc.)** | Bell or glockenspiel — bright, ringing | Single accented notes per metal atom |
| **Other / rare elements** | Synthesizer pad with unique wavetable per element | One-shot textural accent |

**Volume and density rules:**

- Each element's instrument voice has a volume proportional to `count(element) / total_heavy_atoms`.
- Elements appearing only once produce sparse, accented single notes.
- Dominant elements (e.g., carbon in organic molecules) form continuous melodic or harmonic lines.
- The total number of unique instruments in the mix equals the number of distinct elements, directly tying chemical diversity to timbral diversity.

---

### 3. Ring Systems → Bass Loops & Rhythmic Layers

Ring structures are one of the richest sources of musical information. Each ring detected in the molecule generates an independent, looping bass/rhythm layer that repeats throughout the composition.

#### Ring Detection

Rings are identified using RDKit's Smallest Set of Smallest Rings (SSSR). Each ring is analyzed individually and contributes a separate loop layer.

#### Per-Ring Mapping

| Ring Feature | Audio Parameter | Mapping Logic |
|---|---|---|
| **Ring size** (atom count) | **Loop length** (in beats) | 3-membered ring → 3-beat loop, 4 → 4, 5 → 5, 6 → 6, 7+ → 7-beat (odd time feel). Directly encodes ring geometry into rhythmic feel. |
| **Ring atom composition** | **Loop pitch content** | Carbon-only ring → root + fifth (power chord feel). Each heteroatom in the ring adds a chromatic note offset: N → +3 semitones, O → +7, S → +10, etc. More heteroatoms = more dissonant loop. |
| **Aromatic vs. saturated** | **Loop articulation** | Aromatic: sustained, legato, with reverb (smooth, resonant). Saturated: staccato, punchy, dry (percussive). |
| **Ring position in molecule** | **Stereo panning** | Rings are spaced across the stereo field from left to right based on their order in the canonical SMILES traversal. |
| **Substituents on ring** | **Loop modulation** | Number of substituents on the ring controls LFO rate on the loop filter: more substituents → faster filter sweep → more animated loop. |

#### Multi-Ring Behavior

| Multi-Ring Feature | Audio Parameter | Mapping Logic |
|---|---|---|
| **Total ring count** | **Number of simultaneous bass loop layers** | 1 ring = 1 loop. 5 rings = 5 independent loops layered. Each is individually mixed. |
| **Fused rings** (shared edges) | **Polyrhythmic interaction** | Fused rings share a common downbeat but have their individual loop lengths, creating polyrhythmic textures. E.g., a fused 5+6 system creates a 5-against-6 polyrhythm. |
| **Spiro junctions** | **Syncopation accents** | The shared atom in a spiro junction is rendered as a sforzando (accented) hit on both loops simultaneously. |
| **Ring system clusters** | **Bass register grouping** | Isolated ring systems are separated by octave: first system in bass register, second one octave up, etc. |

**Musical rationale:** Benzene produces a smooth, 6-beat legato bass loop. Naphthalene layers two fused 6-beat loops with a polyrhythmic lock. Cyclopropenyl cation gives a tense, 3-beat staccato pattern. A steroid skeleton (four fused rings of different sizes) produces a complex, interlocking rhythmic foundation.

---

### 4. Aromaticity → Drone Pads & Reverb

Aromaticity represents electron delocalization — a smeared, resonant quality — and is mapped to sustained sonic textures.

| Aromatic Feature | Audio Parameter | Mapping Logic |
|---|---|---|
| **Number of aromatic rings** | **Drone pad layers** | Each aromatic ring adds a sustained synth pad layer tuned to the root note of its ring loop (from §3). Zero aromatic rings → no drone. |
| **Total aromatic atom count** | **Reverb depth (wet/dry mix)** | More aromatic atoms → wetter, more reverberant overall mix. 0 aromatic atoms → dry, intimate mix. 20+ → cathedral-like wash. |
| **Aromatic atom fraction** (aromatic / total) | **Pad volume** | Higher fraction → pads dominate the mix, creating an ambient, enveloping feel. Low fraction → pads sit quietly beneath other layers. |
| **Conjugation length** (longest conjugated path) | **Pad sustain / release time** | Longer conjugation → longer sustain and slower release. Short conjugation → quicker pad decay. |
| **Heteroaromatic rings** (pyridine, furan, etc.) | **Pad detuning / dissonance** | Each heteroatom in an aromatic ring detunes the pad by a small interval (±5–15 cents per heteroatom), creating an eerie, shimmering quality. |

**Musical rationale:** Benzene has a single warm pad with moderate reverb. Coronene (7 fused aromatic rings) generates a massive, shimmering drone wall. Aliphatic molecules like cyclohexane are completely dry with no pads.

---

### 5. Main Carbon Chain → Lead Melody

The longest carbon chain (longest path in the molecular graph, considering only C–C bonds) defines the primary melodic line — the "voice" of the molecule.

| Chain Feature | Audio Parameter | Mapping Logic |
|---|---|---|
| **Chain length** (atom count) | **Melody phrase length** | Each carbon in the chain maps to one melodic note. A 2-carbon chain (ethane) = 2 notes. A 20-carbon chain = 20-note phrase. |
| **Carbon index in chain** | **Pitch (scale degree)** | The melody walks through the chosen scale (from §1). Each carbon steps to the next scale degree. Direction (up/down) is influenced by branch points: a branch triggers a direction change. |
| **Chain saturation** | **Legato vs. staccato** | Fully saturated chain segments → legato, connected notes. Unsaturated segments (double/triple bonds within chain) → staccato or accented notes at those positions. |
| **Chain straightness** | **Pitch range (octave span)** | A perfectly linear chain stays within 1 octave. Chains with many direction changes (from branches) span 2–3 octaves. |
| **Position of heteroatoms in chain** | **Chromatic passing tones** | If the longest path includes heteroatoms (e.g., in ethers, amines), those positions insert chromatic notes outside the scale, adding tension. |

**Melody contour algorithm:**

1. Start on the root note.
2. Walk the longest chain atom by atom.
3. At each C: step up one scale degree.
4. At each branch point: reverse direction (up → down or down → up).
5. At each heteroatom: insert chromatic note, then return to scale.
6. At each double bond: add an accent / staccato mark.
7. At each triple bond: hold note for 2x duration.

**Musical rationale:** Methane produces a single sustained note. Butane gives a simple 4-note ascending phrase. Squalene (C₃₀ chain with branches) produces an expansive, winding melody across multiple octaves.

---

### 6. Branches & Ramifications → Counter-Melodies & Harmonies

Every branch (substituent off the main chain) generates a secondary melodic voice that harmonizes with or counterpoints the main melody.

| Branch Feature | Audio Parameter | Mapping Logic |
|---|---|---|
| **Number of branches** | **Number of counter-melody voices** | Each branch = one additional melodic voice. 0 branches → solo melody (monophonic). 5 branches → 5-part harmony. |
| **Branch length** (atoms) | **Counter-melody note count** | Short branches (1–2 atoms): single grace note or trill. Medium (3–5): short motif. Long (6+): full counter-melodic phrase. |
| **Branch position on main chain** | **Counter-melody entry time** | Branches near the start of the chain enter early in the composition. Branches near the end enter late. Creates staggered polyphony. |
| **Branch depth** (nested branches) | **Harmonic interval from root** | Depth 1 (direct substituent): harmony at a third or fifth. Depth 2 (branch off a branch): harmony at a seventh or ninth. Depth 3+: increasingly dissonant intervals (11ths, tritones). |
| **Branch atom composition** | **Counter-melody instrument** | Branch instrument follows the element→instrument mapping from §2, using the most prevalent element in the branch. |
| **Symmetric branches** (identical substituents) | **Unison / octave doubling** | Identical branches produce the same counter-melody in unison or octave doubling, reinforcing that voice. |

**Musical rationale:** Neopentane (C with 4 methyl branches) produces 4 identical grace notes in unison — a single punctuated chord. A dendrimer generates a cascade of staggered, increasingly dissonant counter-melodies. A linear unbranched alkane is a pure solo melody.

---

### 7. Bond Types → Articulation & Dynamics

Bond types modulate *how* notes are played rather than *which* notes are played.

| Bond Feature | Audio Parameter | Mapping Logic |
|---|---|---|
| **Single bond (C–C)** | **Normal articulation** | Default note duration (one beat subdivision). Standard velocity. |
| **Double bond (C=C)** | **Pitch bend / slide** | A pitch bend (portamento) is applied between the two notes connected by a double bond. Creates a "sliding" effect. |
| **Triple bond (C≡C)** | **Sustained note + vibrato** | The note at a triple bond is held for 2× normal duration with added vibrato. Creates a "singing" quality. |
| **Aromatic bond** | **Legato with soft attack** | Notes in aromatic systems have zero attack time and overlap (true legato), sonically reflecting delocalization. |
| **Double bond fraction** (C=C / total bonds) | **Global pitch-bend depth** | Higher fraction → deeper pitch bends across the whole piece. |
| **Triple bond fraction** | **Global vibrato depth** | Higher fraction → more prominent vibrato on all sustained notes. |
| **Bond order sequence along main chain** | **Rhythmic pattern** | The sequence of bond orders (1, 2, 1, 1, 3, 1, …) along the main chain is mapped to a rhythmic pattern: single=eighth note, double=quarter, triple=dotted quarter. Creates a unique rhythm per molecule. |

**Musical rationale:** Ethane (all single bonds) is rhythmically uniform. Acetylene (triple bond) is a long, vibrato-rich sustained note. Polybutadiene (alternating single/double) produces a characteristic swing/shuffle rhythm.

---

### 8. Functional Groups → Motifs & Accents

Recognized functional groups inject pre-designed musical motifs — short, characteristic melodic or rhythmic figures — into the composition at the position where they occur.

| Functional Group | Musical Motif | Description |
|---|---|---|
| **Hydroxyl (–OH)** | Bright upward arpeggio (2 notes) | A quick major-third jump. "Sparkling" quality. |
| **Carbonyl (C=O)** | Bell strike | A single bell tone with long decay. Solitary, resonant. |
| **Carboxylic acid (–COOH)** | Bell strike + upward arpeggio combined | A bell followed by a sparkle — the combination reflects the composite group. |
| **Amine (–NH₂ / –NHR / –NR₂)** | Warm rising chord | A soft sustained chord that swells. Primary → major chord; secondary → minor; tertiary → diminished. |
| **Ester (–COOR)** | Two-note descending figure | A smooth, "sighing" interval. Refined. |
| **Amide (–CONHR)** | Chord + sustained pad swell | Combines the bell of carbonyl with the warmth of the amine. |
| **Ether (–O–)** | Breathy sustained note | A flute-like tone held across the bond. |
| **Thiol (–SH)** | Low growling accent | A distorted bass stab. Aggressive. |
| **Nitro (–NO₂)** | Sharp staccato burst | A rapid double-hit, percussive. Tension. |
| **Phosphate (–PO₄)** | Marimba roll | A rapid tremolo roll on a mallet instrument. |
| **Halide (C–X)** | Metallic ping | A single high-frequency metallic strike, unique per halogen (see §2). |
| **Aldehyde (–CHO)** | Bell strike + grace note | Slightly brighter and lighter than a ketone's bell. |
| **Ketone (C(=O)C)** | Muted bell strike | Deeper, more muffled bell tone than aldehyde. |
| **Nitrile (–C≡N)** | High-pitched sustained whine | A thin, tense, sustained high note with slow vibrato. |
| **Sulfoxide (S=O)** | Distorted pad swell | Low-mid growl that swells. |
| **Alcohol chain (polyol)** | Cascading arpeggios | Multiple hydroxyl motifs staggered in quick succession. |

**Motif placement:** Motifs are inserted at the temporal position corresponding to their location along the main chain traversal. Multiple functional groups produce overlapping motifs.

**Motif stacking:** When a functional group appears multiple times (e.g., glucose with 5 –OH groups), each instance triggers its motif in sequence, creating a characteristic repeating pattern. The repetition count itself becomes a musical feature.

---

### 9. Stereochemistry → Spatial Audio & Phrasing

Stereochemical information encodes 3D spatial arrangement, which maps naturally to spatial audio parameters and melodic direction.

| Stereochemical Feature | Audio Parameter | Mapping Logic |
|---|---|---|
| **Chiral centers (R/S)** | **Stereo panning** | R-configured centers pan their associated notes to the right channel. S-configured centers pan left. A molecule with both R and S centers produces a composition that bounces between left and right. |
| **Number of stereo centers** | **Pan movement frequency** | More stereo centers → more frequent left-right movement in the mix. Creates a dynamic, spatial feel. |
| **E/Z (cis/trans) isomerism** | **Melodic phrase direction** | E (trans) geometry → ascending melodic phrase at that position. Z (cis) → descending phrase. |
| **Meso compounds** (internal symmetry) | **Stereo width collapse** | Meso compounds narrow the stereo field toward center, reflecting their internal compensation of chirality. |
| **Axial chirality** | **Slow auto-pan sweep** | A continuous slow panning sweep rather than discrete L/R jumps. |
| **No stereochemistry** | **Center-panned, static** | Achiral molecules are mixed in mono center. Simple and static spatial image. |

**Musical rationale:** D-glucose and L-glucose produce the same notes but mirrored in the stereo field. A molecule with alternating R/S centers creates a hypnotic ping-pong panning effect. Meso-tartaric acid collapses to center.

---

### 10. Atom Connectivity & Graph Topology → Structure & Form

The molecular graph (atoms as nodes, bonds as edges) defines the high-level musical form — how sections are organized and how complexity unfolds.

| Graph Feature | Audio Parameter | Mapping Logic |
|---|---|---|
| **Graph diameter** (longest shortest path) | **Melodic pitch range (octaves)** | Diameter 1–3 → single octave. 4–8 → two octaves. 9+ → three octaves. Wide-ranging molecules produce wide-ranging melodies. |
| **Average node degree** | **Note overlap / polyphony** | Low average degree (chain-like): monophonic or 2-voice. High average degree (highly connected): 4+ voice polyphony. |
| **Maximum node degree** | **Climax intensity** | The atom with the highest connectivity marks the musical climax point: maximum volume, most voices active, most dissonance. |
| **Number of connected components** | **Number of distinct sections** | Disconnected fragments (salts, complexes) produce separate musical sections with silence between them. E.g., NaCl → two short contrasting phrases. |
| **Graph symmetry** (automorphism group size) | **Repetition and palindrome** | High symmetry → musical phrases are repeated or palindromic (played forward then backward). Low symmetry → through-composed, no repeats. |
| **Bridges** (bonds whose removal disconnects the graph) | **Dramatic pauses** | Bridge bonds in the molecular graph insert brief silences (rests) in the audio, creating phrasing and breath. |
| **Wiener index** (sum of all shortest paths) | **Reverb pre-delay** | Higher Wiener index (more "spread out" molecule) → longer reverb pre-delay, creating a sense of spaciousness. |

**Musical rationale:** Adamantane (highly symmetric cage) produces a palindromic, repeating composition. A long-chain fatty acid (high diameter, low branching) produces a lyrical, wide-ranging solo. A salt pair (disconnected graph) produces a call-and-response form.

---

### 11. Atom Recurrence & Elemental Distribution → Rhythmic Patterns

The frequency and distribution of each element type across the molecule generates rhythmic patterns layered onto the instrument voices from §2.

| Distribution Feature | Audio Parameter | Mapping Logic |
|---|---|---|
| **Element frequency** (count per type) | **Note repetition rate** | An element appearing 10 times generates 10 notes for its instrument voice, evenly distributed across the composition duration. Rare elements get sparse, accented hits; common ones get continuous patterns. |
| **Element clustering** (adjacent same-element atoms) | **Note grouping / tuplets** | Consecutive same-element atoms (e.g., –S–S– disulfide) produce grouped notes: duplets, triplets, etc. A sulfur cluster of 3 → a triplet figure on the S-instrument. |
| **Element periodicity** (regular spacing) | **Rhythmic regularity** | If an element appears at regular intervals in the graph traversal (e.g., every 4th atom), its rhythm is perfectly periodic (on the beat). Irregular spacing → syncopated rhythm. |
| **Element ratio** (e.g., C:O, C:N) | **Rhythmic interplay between instrument pairs** | A 2:1 C:O ratio creates a 2-against-1 polyrhythm between the piano and flute voices. A 3:2 ratio → a characteristic 3:2 hemiola. |
| **Unique element count progression** (traversal order) | **Timbral evolution** | As the molecular graph is traversed and new elements are encountered, their instruments "enter" one by one. A molecule with early diversity (many elements near the start) has all instruments present quickly; one with late diversity builds slowly. |

**Musical rationale:** Glucose (C₆H₁₂O₆) has a C:O ratio of 1:1, producing interlocking equal rhythms. Caffeine, with C, H, N, O all in different ratios, generates a complex layered polyrhythm. A pure hydrocarbon is rhythmically simple (only piano and celesta).

---

### 12. Charges & Electronegativity → Expression & Effects

Electronic properties add expressive nuance — the "emotion" layer of the composition.

| Electronic Feature | Audio Parameter | Mapping Logic |
|---|---|---|
| **Formal positive charge (+)** | **Accent / sforzando** | Each positively charged atom produces a sharp accented note at its position in the traversal. Bright, forward. |
| **Formal negative charge (–)** | **Soft sustained swell** | Negatively charged atoms produce a soft, sustained crescendo-decrescendo at their position. Mellow, receding. |
| **Net molecular charge** | **Overall brightness (EQ)** | Cationic molecules: high-frequency EQ boost (brighter). Anionic: low-frequency boost (darker). Neutral: flat EQ. |
| **Zwitterionic** (both + and –) | **Tension/resolution motif** | A dissonant chord (tension) resolving to consonance, placed at the charged positions. |
| **Electronegativity gradient along chain** | **Pitch direction bias** | If the melody traverses from electropositive to electronegative atoms, the pitch trends upward. The reverse → downward. Encodes the "pull" of electrons as pitch attraction. |
| **Dipole moment estimate** (from charge separation) | **Stereo width** | Large dipole → wide stereo spread. Small/zero dipole → narrow mono image. |
| **Radical electrons** (unpaired) | **Noise / distortion burst** | A brief white noise burst or distortion hit at the position of the radical. Unstable, chaotic. |

---

## Output Formats

### MP3 / WAV (Direct Audio)

```bash
chemsymphony -s "CCO" -f mp3 -o ethanol.mp3
chemsymphony -s "CCO" -f wav -o ethanol.wav
```

All audio layers are synthesized, mixed, and rendered to a single audio file. Default sample rate is 44100 Hz, 16-bit. MP3 uses 192 kbps encoding.

### MIDI + JSON Manifest

```bash
chemsymphony -s "c1ccccc1" -f midi -o benzene/
```

This outputs a directory containing:

```
benzene/
├── manifest.json          # Complete mapping metadata
├── master.mid             # Combined MIDI file (all tracks)
├── tracks/
│   ├── 01_lead_melody.mid     # Main chain melody
│   ├── 02_bass_loop_ring1.mid # Ring bass loops (one per ring)
│   ├── 03_counter_melody_1.mid# Branch counter-melodies
│   ├── 04_drone_pad.mid       # Aromaticity pads
│   ├── 05_percussion.mid      # Halogen/element percussion
│   ├── 06_motifs.mid          # Functional group motifs
│   └── ...
└── README.txt             # Human-readable summary
```

**manifest.json** contains:

```json
{
  "smiles_input": "c1ccccc1",
  "smiles_canonical": "c1ccccc1",
  "molecular_formula": "C6H6",
  "features_extracted": {
    "heavy_atom_count": 6,
    "molecular_weight": 78.11,
    "ring_count": 1,
    "aromatic_ring_count": 1,
    "longest_chain": 6,
    "branch_count": 0,
    "stereo_centers": 0,
    "functional_groups": [],
    "...": "..."
  },
  "audio_parameters": {
    "bpm": 88,
    "key": "D",
    "scale": "pentatonic_major",
    "duration_seconds": 3.8,
    "root_note_midi": 62,
    "...": "..."
  },
  "tracks": [
    {
      "name": "lead_melody",
      "file": "tracks/01_lead_melody.mid",
      "channel": 0,
      "instrument": "acoustic_piano",
      "source_feature": "main_carbon_chain",
      "notes": 6
    },
    {
      "name": "bass_loop_ring1",
      "file": "tracks/02_bass_loop_ring1.mid",
      "channel": 1,
      "instrument": "synth_bass",
      "source_feature": "ring_system_0",
      "loop_beats": 6,
      "aromatic": true
    }
  ]
}
```

The JSON manifest makes ChemSymphony outputs fully inspectable and reproducible, and enables downstream tools (DAWs, visualization engines) to work with individual layers.

---

## Examples

### Simple Molecules

| Molecule | SMILES | Audio Character |
|---|---|---|
| **Methane** | `C` | A single sustained piano note. No bass. No pads. Minimal. |
| **Water** | `O` | A single breathy flute note. Dark scale. Very short. |
| **Ethanol** | `CCO` | 3-note ascending piano melody + one –OH sparkle + one ether breath. Light, cheerful. |
| **Acetic acid** | `CC(=O)O` | 2-note melody with bell strike (C=O) + sparkle (OH). Short, punchy. |

### Medium Molecules

| Molecule | SMILES | Audio Character |
|---|---|---|
| **Benzene** | `c1ccccc1` | 6-note melody + single 6-beat legato bass loop + one warm drone pad. Moderate reverb. Smooth and cyclic. |
| **Aspirin** | `CC(=O)OC1=CC=CC=C1C(=O)O` | Longer melody with branches, ester "sigh" motif, bell strikes, one aromatic bass loop + pad. Medium complexity. |
| **Caffeine** | `CN1C=NC2=C1C(=O)N(C(=O)N2C)C` | Multiple fused ring loops (polyrhythmic), drone pads, multiple bell strikes, N/O instrument voices layered. Rich and driving. |

### Complex Molecules

| Molecule | SMILES | Audio Character |
|---|---|---|
| **Cholesterol** | (long SMILES) | 4 fused ring loops of different sizes (polyrhythmic), long winding melody with many direction changes, methyl branch grace notes, one –OH sparkle, dry mix (no aromaticity). Rhythmically complex, tonally clean. |
| **Chlorophyll** | (long SMILES) | Massive: Mg metal bell, 4+ aromatic ring pads, deep reverb, nitrogen brass voices, long ester chains, dense polyphony. Lush, orchestral, expansive. |
| **Taxol** | (long SMILES) | Extreme: 10+ functional group motifs, multiple rings, heavy branching (5+ counter-melodies), stereo panning from chiral centers, fast tempo, chromatic scale, 3-octave melody range. A full symphonic movement. |

---

## Python API

```python
from chemsymphony import ChemSymphony

# Initialize the engine
cs = ChemSymphony()

# Generate audio bytes
audio_bytes = cs.generate("CCO", format="mp3")

# Save to file
cs.generate_to_file("CCO", output="ethanol.mp3", format="mp3")

# Get MIDI tracks + manifest as dict
result = cs.generate("c1ccccc1", format="midi")
# result.manifest  → dict (the JSON manifest)
# result.master    → bytes (combined MIDI)
# result.tracks    → list of (name, bytes) tuples
result.save("benzene/")

# Extract features only (no audio)
features = cs.extract_features("CC(=O)OC1=CC=CC=C1C(=O)O")
print(features.ring_count)         # 1
print(features.functional_groups)  # ['ester', 'carboxylic_acid']
print(features.longest_chain)      # 4
print(features.audio_parameters)   # {'bpm': 105, 'key': 'A', ...}

# Custom parameters
cs.generate_to_file(
    "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
    output="caffeine.wav",
    format="wav",
    bpm=120,        # Override auto BPM
    key="Am",       # Force A minor
    seed=42         # Reproducible
)
```

---

## Configuration

ChemSymphony can be configured via a `chemsymphony.toml` file placed in the working directory or at `~/.config/chemsymphony/config.toml`.

```toml
[audio]
sample_rate = 44100
bit_depth = 16
mp3_bitrate = 192

[mapping]
# Scale overrides
min_bpm = 70
max_bpm = 180
min_duration = 2.0
max_duration = 120.0

# Volume balance between layers (0.0 – 1.0)
[mapping.mix]
lead_melody = 0.85
bass_loops = 0.70
counter_melodies = 0.55
drone_pads = 0.40
percussion = 0.50
motifs = 0.65
effects = 0.30

[mapping.instruments]
# Override element-to-instrument mapping
# carbon = "electric_guitar"
# oxygen = "violin"
```

---

## Package Structure

```
chemsymphony/
├── pyproject.toml
├── README.md
├── LICENSE
├── chemsymphony/
│   ├── __init__.py
│   ├── __main__.py              # Entry point: python -m chemsymphony
│   ├── cli.py                   # CLI argument parsing (click)
│   ├── core.py                  # ChemSymphony main class
│   ├── canonicalize.py          # SMILES canonicalization & validation
│   ├── features/
│   │   ├── __init__.py
│   │   ├── global_props.py      # §1: MW, atom count, heteroatom ratio
│   │   ├── atoms.py             # §2: Atomic composition analysis
│   │   ├── rings.py             # §3: Ring detection & characterization
│   │   ├── aromaticity.py       # §4: Aromatic system analysis
│   │   ├── chains.py            # §5: Longest chain extraction
│   │   ├── branches.py          # §6: Branch detection & characterization
│   │   ├── bonds.py             # §7: Bond type analysis
│   │   ├── functional_groups.py # §8: Functional group recognition
│   │   ├── stereo.py            # §9: Stereochemistry extraction
│   │   ├── topology.py          # §10: Graph topology analysis
│   │   ├── distribution.py      # §11: Element distribution patterns
│   │   └── electronic.py        # §12: Charges & electronegativity
│   ├── mapping/
│   │   ├── __init__.py
│   │   ├── master.py            # Global parameter mapping
│   │   ├── melody.py            # Lead melody generator
│   │   ├── bass.py              # Ring-based bass loop generator
│   │   ├── harmony.py           # Branch-based counter-melodies
│   │   ├── pads.py              # Aromaticity drone pads
│   │   ├── percussion.py        # Element percussion patterns
│   │   ├── motifs.py            # Functional group motif library
│   │   ├── expression.py        # Stereo, dynamics, effects
│   │   └── form.py              # Graph topology → musical form
│   ├── synthesis/
│   │   ├── __init__.py
│   │   ├── midi_builder.py      # MIDI track assembly
│   │   ├── audio_engine.py      # Direct audio synthesis (scipy)
│   │   ├── mixer.py             # Layer mixing & mastering
│   │   └── effects.py           # Reverb, panning, EQ, distortion
│   ├── export/
│   │   ├── __init__.py
│   │   ├── audio_export.py      # MP3/WAV file writing
│   │   ├── midi_export.py       # MIDI file writing
│   │   └── manifest.py          # JSON manifest generation
│   └── config.py                # Configuration loading
└── tests/
    ├── test_canonicalize.py
    ├── test_features.py
    ├── test_mapping.py
    ├── test_synthesis.py
    ├── test_export.py
    └── test_cli.py
```

---

## Dependencies

| Package | Purpose |
|---|---|
| **rdkit** | SMILES parsing, canonicalization, molecular graph analysis, ring detection, functional group matching, stereochemistry |
| **midiutil** | MIDI file generation |
| **numpy** | Waveform generation, signal processing |
| **scipy** | Audio synthesis (oscillators, filters, envelopes) |
| **pydub** | Audio format conversion and mixing |
| **click** | CLI interface |
| **tomli** | Configuration file parsing |

---

## How Complexity Scales

The mapping system is designed so that molecular complexity directly and naturally produces musical complexity:

| Molecular Property | Simple Molecule (e.g., CH₄) | Complex Molecule (e.g., Taxol) |
|---|---|---|
| Atom count | 1 heavy atom → 1 note, short | 113 heavy atoms → long, dense |
| Tempo | Slow (70 BPM) | Fast (170 BPM) |
| Scale | Pentatonic major (5 notes) | Chromatic (12 notes) |
| Instruments | 1 (piano only) | 8+ (piano, brass, flute, bass, bells, percussion…) |
| Ring loops | 0 bass layers | 6+ polyrhythmic bass layers |
| Drone pads | None (dry mix) | Multiple stacked pads, cathedral reverb |
| Melody | 1 note, 1 octave | 30+ notes, 3 octaves |
| Counter-melodies | 0 (monophonic) | 10+ voices (thick polyphony) |
| Motifs | 0 | 15+ overlapping functional group motifs |
| Stereo | Center, static | Dynamic L/R panning from chiral centers |
| Form | Single phrase | Multi-section, with pauses and climax |

---

## License

MIT License. See [LICENSE](LICENSE) for details.

---

*ChemSymphony: every molecule has a voice.*
