Migrating from RelocaTE2
========================

RelocaTE3 covers the same algorithm as RelocaTE2 but with a modern
Python 3 implementation, no external Perl or blat dependency, and a proper
subcommand CLI instead of generated shell scripts.

Flag mapping
------------

.. list-table::
   :header-rows: 1
   :widths: 35 35 30

   * - RelocaTE2 flag
     - RelocaTE3 flag
     - Notes
   * - ``--fq_dir <dir>``
     - ``--r1 <file> --r2 <file>``
     - Explicit file paths instead of a directory glob
   * - ``--te_fasta``
     - ``--te``
     - Same format (FASTA of TE consensus sequences)
   * - ``--genome_fasta``
     - ``--genome``
     -
   * - ``--reference_ins``
     - ``--repeatmasker``
     - RepeatMasker ``.out`` file; enables Steps 0/6 + FP filter
   * - ``--mismatch``
     - ``--mismatch``
     - Default changed from 2 → 0; use ``--mismatch 2`` to match benchmark
   * - ``--len_cut_match``
     - ``--len-cut-match``
     - Hyphens instead of underscores
   * - ``--len_cut_trim``
     - ``--len-cut-trim``
     -
   * - ``--cpu``
     - ``--threads`` / ``-c``
     -
   * - ``--aligner blat/bwa/bowtie2``
     - *(not available)*
     - minimap2 only; ``-k 11 -w 5`` seeding recovers most blat-only hits
   * - ``--sample_name``
     - ``--sample``
     -
   * - ``characterizer.pl`` (Perl)
     - ``relocaTE3 characterize`` or ``relocaTE3 run --genotype``
     - Pure Python via pysam; same homozygous/heterozygous/somatic logic
   * - ``--config``
     - *(removed)*
     - Tool paths are resolved from ``PATH``; no config file needed

Aligner change
--------------

RelocaTE2 defaulted to **blat** for mapping reads to the TE library, which
finds short and divergent matches that minimap2's default presets miss.
RelocaTE3 compensates by using sensitive minimap2 seeds (``-k 11 -w 5 -N 20
-p 0.5``) for the TE-mapping step and (``-k 13 -w 6``) for genome re-mapping.
On the rice Chr3 2 Mb benchmark this recovers ~89 % of sites vs. RelocaTE2's
~98 %; the gap is in reads with very short (<15 bp) flanking portions that
cannot be placed uniquely regardless of aligner.

Output differences
------------------

- Output files use hyphens in file names where RelocaTE2 used underscores.
- The GFF ``attributes`` field is identical so the same downstream parsers work.
- The ``.txt`` summary is tab-delimited (same schema as RelocaTE2).
- A ``<sample>.fullreads.genome.bam`` is now written alongside the flanking BAM
  and is used for false-junction filtering (not present in RelocaTE2).

Multi-sample runs
-----------------

RelocaTE2 supported a ``--fq_dir`` directory of samples via a multiprocessing
pool driven by generated shell scripts.  RelocaTE3 offers the Python API
:func:`~RelocaTE3.pipeline.run_samples` for the same use case:

.. code-block:: python

   from RelocaTE3.ReadLibrary import ReadLibrary
   from RelocaTE3.pipeline import run_samples

   samples = [
       ReadLibrary(["s1_R1.fq.gz", "s1_R2.fq.gz"], "sample1"),
       ReadLibrary(["s2_R1.fq.gz", "s2_R2.fq.gz"], "sample2"),
   ]
   gffs = run_samples(samples, "TE.fa", "genome.fa", "out",
                      sample_threads=2, step_threads=4)

A Nextflow workflow (Phase 8, planned) will provide the same scatter at cluster
scale.
