RelocaTE3
=========

RelocaTE3 identifies transposable element (TE) insertion polymorphisms from
short-read resequencing data at single-base resolution. It is a pure-Python
reimplementation of `RelocaTE2 <https://github.com/JinfengChen/RelocaTE2>`_,
with no Perl. By default it uses RelocaTE2's aligners: **blat** for the TE
search and **bwa aln** for genome placement, with **samtools**, **bedtools**
and ``pysam`` (and ``bcftools`` for excision mode). Each stage's aligner is
selectable: minimap2, bwa, bwa-mem2, bwa aln, bowtie2 or blat. The pixi
environments and the container image ``ghcr.io/stajichlab/relocate3`` carry
every tool (:doc:`installation`).

Run one sample with ``relocaTE3 run-all``, or many on a laptop, workstation
or SLURM cluster with the Nextflow pipeline in the same repository
(:doc:`pipeline`), which gives every sample the same files as ``run-all``
and adds a cohort genotype matrix and a MultiQC report.

**Status:** release candidate 0.3.0rc1 (see ``CHANGELOG.md``). On the rice
Chr3 2 Mb test data, at the shipped defaults, it recovers 196/200 simulated
insertions at precision 1.000, the same as RelocaTE2;
``tests/acceptance_test.py`` gates this. Validated for paired-end short reads.

.. toctree::
   :maxdepth: 2
   :caption: Contents

   installation
   usage
   pipeline
   output
   migration
   api

Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
