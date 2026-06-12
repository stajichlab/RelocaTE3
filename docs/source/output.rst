Output format
=============

Directory layout
----------------

.. code-block:: text

   <outdir>/
     te_containing/
       <sample>.read_repeat_name.txt     read → TE family assignment
     flanking/
       <sample>.left.flankingReads.fq    trimmed flanking reads (R1)
       <sample>.right.flankingReads.fq   trimmed flanking reads (R2)
     te_portions/
       <sample>.five_prime.fa            5′ TE-matching portions
       <sample>.three_prime.fa           3′ TE-matching portions
     genome_aln/
       <sample>.flanking.genome.bam      flanking reads aligned to genome
       <sample>.fullreads.genome.bam     full (untrimmed) reads (FP filter)
       <sample>.reads.genome.bam         original reads (genotyping, --genotype)
     existingTE.bed                      reference TE annotation (--repeatmasker)
     results/
       <sample>.all_nonref_insert.gff    non-reference insertions (GFF3)
       <sample>.all_nonref_insert.txt    non-reference insertions (tab-delimited)
       <sample>.all_ref_insert.gff       reference/shared insertions (--repeatmasker)
       <sample>.all_ref_insert.txt
       <sample>.all_nonref_insert.characTErized.gff   genotyped (--genotype)
       <sample>.all_nonref_insert.characTErized.txt

GFF3 attribute fields
---------------------

Each non-reference insertion feature has the following GFF3 ``attributes``
column entries, preserving the RelocaTE2 convention:

``Name``
    TE family name (from the TE library FASTA header).

``TSD``
    Target-site duplication sequence (the genomic bases duplicated at the
    insertion point).  ``NA`` when no TSD is detectable.

``Note``
    Insertion category: ``Non-reference``, ``Reference-Only``, or ``Shared``.

``Left_junction_reads``
    Number of junction reads supporting the left (5′) end of the insertion.

``Right_junction_reads``
    Number of junction reads supporting the right (3′) end of the insertion.

``Left_support_reads``
    Number of paired-end reads whose mate is in the TE and whose read maps to
    the left flank.

``Right_support_reads``
    Paired-end support reads on the right flank.

Genotyping outputs
    When ``--genotype`` is enabled, genotyping results are written to
    ``<sample>.all_nonref_insert.characTErized.gff`` and
    ``<sample>.all_nonref_insert.characTErized.txt`` (the non-reference insertion
    GFF is unchanged).

Tab-delimited TXT format
------------------------

The ``.txt`` file has one row per insertion with the same fields as the GFF
``attributes`` column, plus the genomic coordinates, written as plain
tab-separated text for easy downstream processing with ``awk`` or pandas.

read_repeat_name.txt
--------------------

A two-column tab-delimited file mapping each TE-containing read name to the
TE family it matched:

.. code-block:: text

   read_name<TAB>TE_family_name

This table is used by the insertion finder and the reference-insertion caller
to assign TE identity to each junction cluster.
