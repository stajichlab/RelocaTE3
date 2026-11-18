"""Reference access behind one interface (linear FASTA first; graphs later).

Long-read code reads the genome only through :class:`ReferenceSpace`, so a
pangenome graph can later sit behind the same interface.
"""

from __future__ import annotations

from pathlib import Path
from typing import Protocol

import pysam

from RelocaTE3.longread.model import Locus, Reference

# IUPAC complements: R<->Y, K<->M, B<->V, D<->H; N, S and W are their own.
_COMP = str.maketrans(
    "ACGTRYKMSWBVDHNacgtrykmswbvdhn", "TGCAYRMKSWVBHDNtgcayrmkswvbhdn"
)


def reverse_complement(seq: str) -> str:
    """Return the reverse complement of ``seq``, keeping case.

    Args:
        seq: DNA sequence over ``ACGT`` and the IUPAC codes, in either case.

    Returns:
        The reverse complement; ``N``, ``S`` and ``W`` stay themselves.
    """
    return seq.translate(_COMP)[::-1]


class ReferenceSpace(Protocol):
    """The only way long-read code reads reference sequence.

    Attributes:
        reference: The reference this space serves.
    """

    reference: Reference

    def sequences(self) -> list[tuple[str, int]]:
        """Return ``(sequence_id, length)`` pairs in reference order."""
        ...

    def fetch(self, locus: Locus) -> str:
        """Return the uppercase bases of ``locus`` on its strand."""
        ...


class LinearFastaReferenceSpace:
    """A :class:`ReferenceSpace` over one indexed FASTA file.

    Attributes:
        reference: The reference this space serves.
    """

    def __init__(self, reference: Reference) -> None:
        """Open the FASTA, creating its ``.fai`` index if missing.

        Args:
            reference: Reference whose ``fasta`` path is opened.

        Raises:
            FileNotFoundError: If the FASTA does not exist.
        """
        path = Path(reference.fasta)
        if not path.is_file():
            raise FileNotFoundError(f"reference FASTA not found: {path}")
        if not Path(str(path) + ".fai").is_file():
            pysam.faidx(str(path))
        self.reference = reference
        self._fasta = pysam.FastaFile(str(path))

    def sequences(self) -> list[tuple[str, int]]:
        """Return ``(sequence_id, length)`` pairs in FASTA file order."""
        return list(zip(self._fasta.references, self._fasta.lengths))

    def fetch(self, locus: Locus) -> str:
        """Return the uppercase bases of ``locus``, reverse-complemented on ``-``.

        Args:
            locus: Interval to fetch; must be on this reference.

        Returns:
            The bases, uppercased (pysam keeps soft-mask case as stored).

        Raises:
            ValueError: If ``locus.reference_id`` is not this reference.
            KeyError: If ``locus.sequence_id`` is not in the FASTA.
        """
        if locus.reference_id != self.reference.id:
            raise ValueError(
                f"locus is on {locus.reference_id!r}, not {self.reference.id!r}"
            )
        if locus.sequence_id not in self._fasta.references:
            raise KeyError(locus.sequence_id)
        seq = self._fasta.fetch(locus.sequence_id, locus.start, locus.end).upper()
        return reverse_complement(seq) if locus.strand == "-" else seq
