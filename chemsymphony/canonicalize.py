"""SMILES canonicalization and validation via RDKit."""

from __future__ import annotations

from dataclasses import dataclass

from rdkit import Chem


class InvalidSmilesError(ValueError):
    """Raised when a SMILES string cannot be parsed."""


@dataclass(frozen=True)
class CanonicalResult:
    """Result of SMILES canonicalization."""

    input_smiles: str
    canonical_smiles: str
    mol: Chem.Mol


def canonicalize(smiles: str) -> CanonicalResult:
    """Validate and canonicalize a SMILES string.

    Returns a ``CanonicalResult`` containing the canonical SMILES and the
    RDKit ``Mol`` object ready for downstream analysis.

    Raises ``InvalidSmilesError`` if the SMILES string is invalid.
    """
    if not smiles or not smiles.strip():
        raise InvalidSmilesError("Empty SMILES string")

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise InvalidSmilesError(f"Invalid SMILES: {smiles!r}")

    canonical = Chem.MolToSmiles(mol, canonical=True)

    # Re-parse the canonical SMILES so the Mol is in canonical form
    mol = Chem.MolFromSmiles(canonical)
    assert mol is not None

    return CanonicalResult(
        input_smiles=smiles,
        canonical_smiles=canonical,
        mol=mol,
    )
