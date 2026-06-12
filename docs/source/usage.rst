Usage
=====

Quick start
-----------

Run the full pipeline for one paired-end sample:

.. code-block:: bash

   relocaTE3 run \
     --r1 reads_1.fq.gz --r2 reads_2.fq.gz \
     --te RiceTE.fa \
     --genome reference.fa \
     --outdir HEG4_out --sample HEG4 \
     --threads 8 --mismatch 2

Add ``--repeatmasker reference.fa.RepeatMasker.out`` to enable reference/shared
insertion calls and false-positive filtering.  Add ``--genotype`` to classify
each insertion as homozygous, heterozygous, or somatic.

Subcommands
-----------

Each pipeline step is also a standalone subcommand so a workflow engine
(e.g. Nextflow) can scatter them across chromosomes or samples.

.. code-block:: bash

   # Step 3 — trim TE-matching portions from reads
   relocaTE3 trim \
     --r1 reads_1.fq.gz --r2 reads_2.fq.gz \
     --te RiceTE.fa \
     --outdir sample_out --sample HEG4

   # Step 4 — re-align flanking reads to the genome
   relocaTE3 align-genome \
     --r1 reads_1.fq.gz --r2 reads_2.fq.gz \
     --genome reference.fa \
     --outdir sample_out --sample HEG4

   # Step 5 — call non-reference insertions
   relocaTE3 find-insertions \
     --bam sample_out/genome_aln/HEG4.genome.bam \
     --genome reference.fa \
     --read-repeat sample_out/te_containing/HEG4.read_repeat_name.txt \
     --outdir sample_out --sample HEG4

   # Steps 0/6 — reference and shared insertions
   relocaTE3 find-reference \
     --bam sample_out/genome_aln/HEG4.genome.bam \
     --repeatmasker reference.fa.RepeatMasker.out \
     --read-repeat sample_out/te_containing/HEG4.read_repeat_name.txt \
     --outdir sample_out --sample HEG4

   # Step 7 — genotype insertions
   relocaTE3 characterize \
     --insertions sample_out/results/HEG4.all_nonref_insert.gff \
     --reads-bam sample_out/genome_aln/HEG4.reads.genome.bam \
     --outdir sample_out --sample HEG4

Common parameters
-----------------

``--mismatch`` *N*
    Mismatches allowed when aligning reads to the TE library (default 0).
    Use 2 to match the sensitivity of RelocaTE2's benchmark settings.

``--len-cut-match`` *N*
    Minimum number of bases a read must match the TE consensus (default 10).

``--len-cut-trim`` *N*
    Minimum length of the trimmed flanking portion to retain a read (default 10).

``--threads`` / ``-c`` *N*
    CPU threads passed to minimap2 (default 1).

``--repeatmasker`` *FILE*
    RepeatMasker ``.out`` for the reference genome.  Enables reference/shared
    insertion calls and drops non-reference calls that overlap a same-family
    reference TE.

``--genotype``
    Map the original reads to the genome and classify each insertion as
    homozygous, heterozygous, or somatic (replaces ``characterizer.pl``).

Python API
----------

The same steps are importable as a library:

.. code-block:: python

   from RelocaTE3.ReadLibrary import ReadLibrary
   from RelocaTE3.pipeline import run_sample, run_samples

   reads = ReadLibrary(["reads_1.fq.gz", "reads_2.fq.gz"], "HEG4")

   # single sample
   gff = run_sample(reads, "RiceTE.fa", "reference.fa", "HEG4_out",
                    mismatch_allowance=2, threads=8)

   # multiple samples in parallel (2 samples at a time, 4 threads each)
   all_reads = [reads, ReadLibrary(["other_1.fq.gz", "other_2.fq.gz"], "HEG5")]
   gffs = run_samples(all_reads, "RiceTE.fa", "reference.fa", "out",
                      sample_threads=2, step_threads=4)
