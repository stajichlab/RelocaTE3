"""Long-read TE insertion detection (milestone L2 onward).

This subpackage aligns long reads to one or more references, turns the
alignments into format-independent placement records, extracts insertion
signatures (CIGAR insertions, split reads, clipped ends) and clusters them
into candidate insertion sites. Every position carries the reference it
belongs to, so nothing here assumes a single reference genome.
"""
