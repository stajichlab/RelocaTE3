"""ReadLibrary for paired-end or single-end read file set."""

from __future__ import annotations

from collections import defaultdict


class ReadLibrary:
    """Represent sequence library typically a single or paired-end FASTQ read files."""

    is_paired = False  # Paired library status
    file_set = [None, None]  # set of files usually just 1 or 2 to store
    name = ""  # Name of the individual strain (eg ReadGroup)

    def __init__(self, fileset: list[str], ind_name: str):
        """Initialize the ReadLibrary."""
        if fileset is not None:
            if len(fileset) > 2:
                raise ValueError("Fileset needs to be either one or two files provided")
            self.is_paired = len(fileset) == 2
            self.file_set = fileset
            self.name = ind_name

    def left(self) -> str:
        """Shortcut to get left read pair."""
        return self.file_set[0]

    def right(self) -> str:
        """Shortcut to get right read pair."""
        return self.file_set[1]


class TrimmedReadLibrary:
    """Read library files processed for Transposon trimmed sequence."""

    def __init__(self, fileset: list[str] = [], ind_name: str = ""):
        """Initialize the ReadLibrary."""
        if fileset is not None:
            if len(fileset) > 2:
                raise ValueError("Fileset needs to be either one or two files provided")
            self.is_paired = len(fileset) == 2
            self.file_set = fileset
            self.name = ind_name
        self.flanking_reads = [
            dict(),
            dict(),
        ]  # reads trimmed of TE sequence and reads matched to the middle of TE
        self.containing_reads = [dict(), dict()]  # entire reads matches to TE
        self.trimmed_coordinates = defaultdict(lambda: defaultdict(lambda: str))
        # TE proportion of reads that matched to TE
        self.read_TE_match = {"five_prime": [None, None], "three_prime": [None, None]}
        self.num_matched_reads = 0
        self.num_flank_reads = 0
        # reads have two or more matches on TE and on same strand.
        # Tandem insertions and Highly false positive.
        self.potential_tandemInserts = [None, None]

    def sync_reads(self, outdir: str):
        """Write processed reads info to file for storage."""

    def _add_flank_read(self, left: bool, read_name: str, start: int, end: int):
        """Add a read which overlaps a TE with coordinates.

        Args:
            left (bool): _description_
            read_name (str): _description_
            start (int): _description_
            end (int): _description_
        """
        lib_index = 1 if not left else 0
        self.flanking_reads[lib_index][read_name] = [start, end]
        if read_name not in self.flanking_reads[1 - lib_index]:
            self.flanking_reads[1 - lib_index][read_name] = [-1, -1]

    class TrimmedRead:
        """Trimmed Read info."""

        def __init__(self):
            """Init library."""
            self.read_name = ""  # Name of the read
            self.target_name = ""  # Name of the TE matched
            self.read_length = 0  # Length of the read (Query)
            self.target_length = 0  # Length of the TE (Target)
            self.start = 0  # Start position of alignment in the read
            self.end = 0  # End position of alignment in the read
            self.target_start = 0  # Start position of alignment in the target
            self.target_end = 0  # End position of alignment in the target
            self.mismatch = 0  # Number of mismatches
            self.strand = ""  # Strand of alignment (in read)
            self.boundary = ""  # Boundary of alignment (in read)
            self.score = 0  # Alignment score
            self.cigar = ""  # CIGAR string
            self.mapq = 0  # Mapping quality
