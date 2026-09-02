Output format
=============

Directory layout
----------------

.. code-block:: text

   <outdir>/
     <sample>.left.bam / <sample>.right.bam    reads aligned to the TE library (map)
     <sample>.repeat.minimap.sorted.bam        flanking reads aligned to genome (align-genome)
     te_containing/
       <sample>.read_repeat_name.txt           read → TE family assignment
       <sample>.left.ContainingReads.fq
       <sample>.right.ContainingReads.fq
     flanking/
       <sample>.left.flankingReads.fq          trimmed flanking reads (5′)
       <sample>.right.flankingReads.fq         trimmed flanking reads (3′)
     te_portions/
       <sample>.five_prime.fa                  5′ TE-matching portions
       <sample>.three_prime.fa                 3′ TE-matching portions
     results/
       <target>.<te_name>.all_nonref_insert.txt                 non-reference insertions (find-insertions)
       <target>.<te_name>.all_nonref_insert.characTErized.gff   genotyped insertions (characterize)
       <target>.<te_name>.all_nonref_insert.characTErized.txt

``<target>`` and ``<te_name>`` come from ``find-insertions --target`` / ``--te-name``.
The ``annotate-ref`` subcommand separately writes ``existingTE.bed`` (reference TE
copies) when run.

GFF3 attribute fields
---------------------

Each non-reference insertion feature has the following GFF3 ``attributes``
column entries, preserving the RelocaTE2 convention:

``Name``
    Single primary TE family name (from the TE library FASTA header).

``TE_family_support``
    Comma-separated ``family=count`` junction-read votes, ordered by decreasing
    support and then family name.

``TE_family_confidence``
    Fraction of informative family votes supporting ``Name``. This measures
    assignment agreement, not insertion-call confidence.

``TE_family_status``
    ``unique`` when only one family received votes, ``dominant`` when the
    primary family has an absolute majority, ``ambiguous`` when no family has
    a majority, or ``unassigned`` when no informative family vote is available.

``TE_supporting_family_support``
    Comma-separated ``family=count`` votes from bracketing supporting mates.
    These are reported separately because a mate links the TE family to the
    locus indirectly and cannot change the junction-driven ``Name``.

``TE_supporting_family_confidence``
    Fraction of informative supporting-mate votes assigned to their most
    supported family.

``TE_supporting_family_status``
    ``unique``, ``dominant``, ``ambiguous``, or ``unassigned``, using the same
    vote rules as the junction-family status.

``TE_family_concordance``
    Relationship between the junction and supporting primary families:
    ``concordant``, ``discordant``, ``junction_only``, ``supporting_only``, or
    ``unassigned``.

``TSD``
    Target-site duplication sequence (the genomic bases duplicated at the
    insertion point).  ``NA`` when no TSD is detectable.

``Note``
    Insertion category. The installed staged pipeline emits ``Non-reference``;
    ``Reference-Only`` / ``Shared`` come from reference/shared calling (planned
    ``find-reference`` subcommand).

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
    The ``characterize`` subcommand writes genotyping results to
    ``<target>.<te_name>.all_nonref_insert.characTErized.gff`` and
    ``<target>.<te_name>.all_nonref_insert.characTErized.txt`` (the input
    non-reference table is unchanged).

Tab-delimited TXT format
------------------------

The ``.txt`` file has one row per insertion with the same fields as the GFF
``attributes`` column, plus the genomic coordinates, written as plain
tab-separated text for easy downstream processing with ``awk`` or pandas. The
family-evidence columns are appended so the established RelocaTE2-compatible
columns retain their positions.

read_repeat_name.txt
--------------------

A two-column tab-delimited file mapping each TE-containing read name to the
TE family it matched:

.. code-block:: text

   read_name<TAB>TE_family_name

This table is used by the insertion finder and the reference-insertion caller
to assign TE identity to each junction cluster.
