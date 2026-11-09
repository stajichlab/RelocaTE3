Installation
============

RelocaTE3 is a Python package (Python 3.10 or newer) that shells out to
external aligners and samtools. There are three ways to get it with its tools:
pixi environments from a clone, the container image, or pip into an
environment that already has the tools.

External tools
--------------

At the shipped defaults (``--te-aligner blat``, ``--genome-aligner bwaaln``,
matching RelocaTE2), ``relocaTE3 run-all`` calls:

- `blat <https://genome.ucsc.edu/FAQ/FAQblat.html>`_: the default TE aligner.
  BLAT is free for academic, nonprofit and personal use only; commercial users
  need a license from its authors.
- `bwa <https://github.com/lh3/bwa>`_: the default genome aligner (``bwa aln``).
- `samtools <https://www.htslib.org/>`_ ≥ 1.17.
- `minimap2 <https://github.com/lh3/minimap2>`_ ≥ 2.24: ``index-genome`` builds
  its index, and ``--te-aligner`` / ``--genome-aligner minimap2`` use it.
- `seqtk <https://github.com/lh3/seqtk>`_: FASTQ to FASTA for BLAT. Without
  it RelocaTE3 converts in Python, which is much slower on a large library.
- `bcftools <https://samtools.github.io/bcftools/>`_: excision mode only
  (``characterize -x``, ``--excision``).

The other aligners (``bwa-mem2``, ``bowtie2``) are needed only when you
select them. The pixi environments and the container image include all of the
above, plus ``bedtools``.

With pixi, from a clone
-----------------------

`pixi <https://pixi.sh>`_ installs pinned versions of every tool from
``pixi.lock``. ``pixi.toml`` defines four environments:

.. list-table::
   :header-rows: 1
   :widths: 15 45 40

   * - Environment
     - Contents
     - Use it for
   * - ``blat``
     - RelocaTE3 and the core tools plus ``blat`` (linux-64 only)
     - ``relocaTE3 run-all`` with the default aligners
   * - ``pipeline``
     - ``blat`` plus MultiQC (linux-64 only)
     - The Nextflow pipeline's ``pixi`` profile; the container image is built
       from it
   * - ``nf``
     - Nextflow ≥ 26.04 and nf-test only
     - Running the Nextflow pipeline and its tests
   * - ``default``
     - RelocaTE3, the core tools and the plotting stack; **no blat**
     - The test suite (``pixi run test``), the docs, the real-rice validation
       reports

``blat`` is kept out of ``default`` because its bioconda build pins
``zlib 1.2.11``, which the plotting stack will not solve against.

.. code-block:: bash

   git clone https://github.com/stajichlab/RelocaTE3
   cd RelocaTE3
   pixi install -e blat
   pixi run -e blat relocaTE3 --help
   pixi run -e blat relocaTE3 run-all ...

In ``default``, ``relocaTE3`` runs only with a non-default TE aligner (for
example ``--te-aligner minimap2``), unless a ``blat`` is on ``PATH``. The
other aligners find far fewer insertions than ``blat``.

Container image
---------------

``ghcr.io/stajichlab/relocate3`` holds RelocaTE3 and the exact tool versions
of the ``pipeline`` environment. The Nextflow pipeline uses it with
``-profile apptainer``, ``singularity`` or ``docker``
(:doc:`pipeline`). Pin a fixed tag: a release version (``0.3.0rc1``) or
``sha-<short>``; ``dev`` and ``pr-<n>`` move with every push.

.. code-block:: bash

   apptainer exec docker://ghcr.io/stajichlab/relocate3:0.3.0rc1 relocaTE3 --help

Nextflow pipeline
-----------------

The pipeline (``main.nf``) needs Nextflow 26.04 or newer and one source of
tools: the container image or, from a clone, the ``pipeline`` environment
(``-profile pixi``). Get Nextflow from the ``nf`` environment
(``pixi run -e nf nextflow ...``) or a module (on UCR HPCC:
``module load nextflow/26.04.6 apptainer/1.5.3``). See :doc:`pipeline`.

With pip
--------

pip installs RelocaTE3 and its Python dependency (``pysam``) only. Install
into an environment that already has ``blat``, ``bwa``, ``samtools``,
``minimap2``, ``seqtk`` and, for excision mode, ``bcftools`` on ``PATH``:

.. code-block:: bash

   pip install .        # or: pip install -e . for development
   relocaTE3 --help

On an HPC cluster, load the tools as modules first, if your cluster provides
them under these names (check with ``module avail``):

.. code-block:: bash

   module load blat bwa samtools minimap2 seqtk bcftools
   pip install --user .
