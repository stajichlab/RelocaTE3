Usage
=====

Quick start
-----------

RelocaTE3 runs as a sequence of staged subcommands.  The ``run`` subcommand
performs **TE-read identification and flank generation only** (map + trim,
steps 2-3) — it is *not* the complete pipeline.  A full single-sample analysis
chains the stages:

.. code-block:: bash

   # 1. index the reference genome (one-time)
   relocaTE3 index-genome -g reference.fa

   # 2-3. identify TE-containing reads and emit flanking reads (map + trim)
   relocaTE3 run \
     --left reads_1.fq.gz --right reads_2.fq.gz \
     --te-library RiceTE.fa \
     --name HEG4 --outdir HEG4_out \
     --threads 8 --mismatch 2

   # 4. re-align trimmed flanking reads to the genome
   relocaTE3 align-genome \
     -g reference.fa \
     -f HEG4_out/flanking/HEG4.left.flankingReads.fq HEG4_out/flanking/HEG4.right.flankingReads.fq \
     -1 reads_1.fq.gz -2 reads_2.fq.gz \
     -n HEG4 -o HEG4_out --threads 8

   # 5. call non-reference insertions
   relocaTE3 find-insertions \
     -b HEG4_out/HEG4.repeat.minimap.sorted.bam \
     --read-repeat HEG4_out/te_containing/HEG4.read_repeat_name.txt \
     --tsd TTA --target ALL --name HEG4 --outdir HEG4_out --te-name mping \
     --fullreads-bam HEG4_out/genome_aln/HEG4.fullreads.genome.bam \
     --reference-ins reference.fa.RepeatMasker.out

   # 7. genotype the insertions from a reads-to-genome BAM/CRAM
   relocaTE3 characterize \
     -s HEG4_out/results/ALL.mping.all_nonref_insert.txt \
     -b reads_to_genome.bam -g reference.fa \
     -o HEG4_out/results

Subcommands
-----------

Each pipeline step is a standalone subcommand so a workflow engine (e.g.
Nextflow) can scatter them across chromosomes or samples.  The installed
subcommands are ``index-genome``, ``map``, ``trim``, ``run`` (map + trim),
``annotate-ref``, ``align-genome``, ``find-insertions``, and ``characterize``.
Run ``relocaTE3 <command> --help`` for each command's flags.

.. code-block:: bash

   # Step 2 — align reads to the TE library (separate map + trim)
   relocaTE3 map \
     -l reads_1.fq.gz -r reads_2.fq.gz \
     -T RiceTE.fa -n HEG4 -o HEG4_out

   # Step 3 — trim TE-matching portions from the TE-library BAMs
   relocaTE3 trim \
     -b HEG4_out/HEG4.left.bam HEG4_out/HEG4.right.bam \
     -n HEG4 -o HEG4_out

   # CRAM input to characterize requires the reference FASTA
   relocaTE3 characterize \
     -s HEG4_out/results/ALL.mping.all_nonref_insert.txt \
     -b original_reads.cram -g reference.fa \
     -o HEG4_out/results

Common parameters
-----------------

``--mismatch`` *N*
    Mismatches allowed when aligning reads to the TE library, and when
    comparing reads to the genome (default 2, matching RelocaTE2's
    ``--mismatch`` and ``--mismatch_junction``).

``--min-match`` *N*
    Minimum number of bases a read must match the TE consensus
    (``run`` / ``trim``; default 10).

``--min-trimmed`` *N*
    Minimum length of the trimmed flanking portion to retain a read
    (``run`` / ``trim``; default 10).

``--threads`` *N*
    CPU threads passed to the aligner (default 1).

``--tsd`` *MOTIF*
    (``find-insertions``) Target-site duplication motif (e.g. ``TTA``), a fixed
    wildcard (e.g. ``...``), or ``UNK`` to infer each TSD (default ``UNK``,
    which is what RelocaTE2 always uses).

``--reference-ins`` *FILE*
    (``find-insertions``) RepeatMasker ``.out`` or BED of reference TE copies.
    Drops non-reference calls that overlap a same-family reference TE.

``-x`` / ``--excision``
    (``characterize``) Also search for excision events that leave a footprint.

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
