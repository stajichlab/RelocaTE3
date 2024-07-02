"""ReadLibrary for paired-end or single-end read file set."""

from __future__ import annotations


class ReadLibrary():
    """Represent sequence library typically a single or paired-end FASTQ read files."""

    is_paired = False                   # Paired library status
    file_set = [None, None]              # set of files usually just 1 or 2 to store
    name = ""                           # Name of the individual strain (eg ReadGroup)

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


class TrimmedReadLibrary():
    """Read library files processed for Transposon trimmed sequence."""

    name = ""            # Name of the individual strain (eg ReadGroup)
    is_paired = False    # Paired library status
    flanking_reads = [dict(), dict()]     # reads trimmed of TE sequence and reads matched to the middle of TE
    containing_reads = [dict(), dict()]   # all reads have matches to TE
    read_TE_match = {'five_prime': [None, None],   # TE proportion of reads that matched to TE
                     'three_prime': [None, None]}
    num_matched_reads = 0
    num_flank_reads = 0

    # reads have two or more matches on TE and on same strand.
    # Tandem insertions and Highly false positive.
    potential_tandemInserts = [None, None]

    def __init__(self, ind_name: str):
        """Initialize the TrimmedReadLibrary."""
        self.name = ind_name

    def sync_reads(self):
        """Write processed reads info to file for storage."""
        print('writing')

    def _add_flank_read(self, left: bool, read_name: str, start: int, end: int):
        lib_index = 0
        if not left:
            lib_index = 1

        self.flanking_reads[lib_index][read_name] = [start, end]
        if read_name not in self.flanking_reads[1 - lib_index]:
            self.flanking_reads[1 - lib_index][read_name] = [-1, -1]
