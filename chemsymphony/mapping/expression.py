"""§9+12 mapping: Stereochemistry + charges → spatial audio & expression."""

from __future__ import annotations

from dataclasses import dataclass, field

from chemsymphony.config import Config
from chemsymphony.features import MolecularFeatures


@dataclass
class ExpressionParams:
    """Expression and effects parameters applied globally or per-note."""

    # Stereo
    stereo_events: list[dict] = field(default_factory=list)
    stereo_width: float = 1.0
    is_meso: bool = False

    # EQ / brightness
    eq_brightness: float = 0.0  # -1 = dark (anionic), +1 = bright (cationic)

    # Charge accents
    charge_accents: list[dict] = field(default_factory=list)

    # Radicals → noise bursts
    noise_bursts: list[int] = field(default_factory=list)

    # Reverb
    reverb_depth: float = 0.0
    reverb_predelay: float = 0.0


def generate_expression(feat: MolecularFeatures, cfg: Config) -> ExpressionParams:
    """Compute expression parameters from stereochemistry and electronic features."""
    expr = ExpressionParams()

    # §9: Stereo panning from chiral centers
    for center in feat.chiral_centers:
        pan = 0.7 if center["config"] == "R" else -0.7 if center["config"] == "S" else 0.0
        expr.stereo_events.append({
            "atom_idx": center["atom_idx"],
            "pan": pan,
        })

    # Meso → narrow stereo
    expr.is_meso = feat.is_meso
    if feat.is_meso:
        expr.stereo_width = 0.2
    elif feat.stereo_center_count > 0:
        expr.stereo_width = 1.0
    else:
        expr.stereo_width = 0.5

    # §12: Net charge → EQ brightness
    if feat.net_charge > 0:
        expr.eq_brightness = min(1.0, feat.net_charge * 0.3)
    elif feat.net_charge < 0:
        expr.eq_brightness = max(-1.0, feat.net_charge * 0.3)

    # Charge accents
    for charge in feat.formal_charges:
        expr.charge_accents.append({
            "atom_idx": charge["atom_idx"],
            "charge": charge["charge"],
            "accent": "sforzando" if charge["charge"] > 0 else "swell",
        })

    # Radicals → noise positions
    expr.noise_bursts = list(feat.radical_electrons)

    # Reverb from aromaticity
    expr.reverb_depth = min(1.0, feat.aromatic_atom_count / 20.0)

    # Reverb pre-delay from Wiener index
    expr.reverb_predelay = min(100.0, feat.wiener_index / 100.0)

    return expr
