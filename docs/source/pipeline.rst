Nextflow pipeline
=================

The repository is also a Nextflow pipeline (nf-core layout). It runs
RelocaTE3 on many samples at once and adds a cohort table and a summary
report.

Each task runs one ``relocaTE3`` subcommand with the settings
``relocaTE3 run-all`` uses. So every sample's files are identical to
``run-all --genotype`` with the same settings. (A sample merged from several
lanes has no ``run-all`` equivalent.) The workflow holds no detection logic.
``tests/nextflow/parity_chr3.sh`` checks the identity on the Chr3 test data.

What the pipeline adds over ``run-all``:

- samples run in parallel, and the genome is indexed once;
- the genotyping alignment runs alongside the TE search;
- the BLAT TE search of one sample is split over many tasks
  (`BLAT scatter`_);
- ``cohort/`` tables and one MultiQC report over all samples;
- on UCR HPCC, each task is a SLURM job.

.. contents:: On this page
   :local:
   :depth: 1


Quick start
-----------

Run the bundled Chr3 test data (2 Mb of rice chromosome 3, 200 simulated
insertions) straight from GitHub with Apptainer:

.. code-block:: bash

   nextflow run stajichlab/RelocaTE3 -r v0.3.0rc1 -profile test,apptainer \
     --container_tag 0.3.0rc1

Results go to ``results/nextflow_test/`` in the launch directory. The run
takes a few minutes after the image is pulled.

.. note::

   For the development version, use ``-r main --container_tag dev``. The
   ``dev`` image is built from ``main``:

   .. code-block:: bash

      nextflow run stajichlab/RelocaTE3 -r main -profile test,apptainer \
        --container_tag dev

Your own data:

.. code-block:: bash

   nextflow run stajichlab/RelocaTE3 -r v0.3.0rc1 -profile apptainer \
     --container_tag 0.3.0rc1 \
     --input samples.csv \
     --genome reference.fa \
     --te_library mping.fa --te_name mping \
     --repeatmasker reference.fa.RepeatMasker.out \
     --outdir results

``--input``, ``--genome``, ``--te_library`` and ``--outdir`` are required.
``--repeatmasker`` is optional.


Requirements
------------

- **Nextflow 26.04 or newer.** The pipeline uses the nf-schema 3.0.0 plugin,
  which needs 26.04. Older Nextflow stops at start.
- **One way to get the tools**, chosen with ``-profile``:

  - Apptainer or Singularity (recommended on clusters);
  - Docker;
  - pixi, from a clone of the repository (for development).

The container image ``ghcr.io/stajichlab/relocate3`` holds the exact tool
versions RelocaTE3 is validated with, plus MultiQC. It is built from the
``pipeline`` environment of ``pixi.lock``. The download is about 340 MB, about
1.1 GB unpacked.

On UCR HPCC both are modules:

.. code-block:: bash

   module load nextflow/26.04.6 apptainer/1.5.3

From a clone, without containers (for development):

.. code-block:: bash

   git clone https://github.com/stajichlab/RelocaTE3 && cd RelocaTE3
   pixi install -e pipeline -e nf      # task tools + Nextflow and nf-test
   pixi run -e nf nextflow run . -profile test,pixi     # Chr3 test data, about 2 min
   pixi run -e nf nextflow run . -profile pixi \
     --input samples.csv --genome reference.fa --te_library RiceTE.fa \
     --repeatmasker reference.fa.RepeatMasker.out --te_name mping --outdir results

The ``pipeline`` environment is the ``blat`` environment plus MultiQC, with
identical tool builds; the image is built from it. On the Chr3 test data the
pipeline's files are identical to ``run-all``'s with pixi, Docker and
Apptainer (``tests/nextflow/parity_chr3.sh``).


Sample sheet
------------

``--input`` is a CSV file with a header line:

.. code-block:: text

   sample,fastq_1,fastq_2,genotype_bam
   HEG4,/data/HEG4_L1_R1.fq.gz,/data/HEG4_L1_R2.fq.gz,
   HEG4,/data/HEG4_L2_R1.fq.gz,/data/HEG4_L2_R2.fq.gz,
   A123,/data/A123_R1.fq.gz,/data/A123_R2.fq.gz,/data/A123.cram
   S7,/data/S7.fq.gz,,

Columns:

- ``sample`` (required): the sample name. It becomes the output directory
  ``<outdir>/<sample>/`` and part of file names.
- ``fastq_1`` (required): read 1 FASTQ (``.fq``, ``.fastq``, optionally
  ``.gz``).
- ``fastq_2`` (optional): read 2. Leave it empty for single-end reads.
- ``genotype_bam`` (optional): an existing BAM or CRAM of this sample's reads
  on the same genome. It is used for genotyping instead of aligning the reads
  again. Its index must sit beside it: ``<file>.bai`` or ``<file>.csi`` for a
  BAM, ``<file>.crai`` for a CRAM (``samtools index``).

Rules:

- **Lanes.** Rows with the same ``sample`` are lanes of one sample. They are
  concatenated in sheet order. All lanes of a sample must be either paired or
  single-end, and must give the same ``genotype_bam`` (or none).
- **Names** may contain only letters, digits, ``.``, ``_`` and ``-``, and must
  start with a letter or digit.
- **Reserved names:** ``cohort``, ``multiqc`` and ``pipeline_info``. These are
  output directories beside the samples' directories.
- **Paths.** Relative paths are resolved against the directory Nextflow is
  launched from, not the sheet's directory. Use absolute paths.
- A file name may not contain spaces.

The sheet is checked before any task runs. Missing files, bad names and
missing indexes stop the run with a message naming the sample.

A sample with a ``genotype_bam`` is always genotyped, even with
``--genotype false``. This is what ``run-all --genotype-bam`` does.


Parameters
----------

Settings use ``run-all``'s names in snake case, with the same defaults, so the
defaults reproduce RelocaTE2. The exception is ``--genotype``, which is on by
default in the pipeline (``run-all`` needs ``--genotype``).

Pass parameters on the command line (``--size 300``) or in a YAML or JSON file
with ``-params-file`` (see `For agents and scripts`_). On/off settings take
``true`` or ``false``: ``--genotype false`` works. A bare flag
(``--excision``) means ``true``. Any other word stops the run.

``nextflow run stajichlab/RelocaTE3 --help`` prints the parameters.
``--help <parameter>`` describes one in detail.

The tables below are generated from ``nextflow_schema.json``.

.. BEGIN GENERATED PARAMETERS (docs/source/_gen_params.py from nextflow_schema.json; do not edit by hand)

Input/output options
~~~~~~~~~~~~~~~~~~~~

Sample sheet, shared reference inputs and where results are written.

The sample sheet lists one row per sample (or per lane: rows sharing a sample name are merged). The genome, TE library and RepeatMasker annotation are shared by all samples.

.. list-table::
   :header-rows: 1
   :widths: 22 10 16 52

   * - Parameter
     - Type
     - Default
     - Description
   * - ``--input``
     - string
     - (none)
     - Sample sheet (CSV): sample,fastq_1,fastq_2[,genotype_bam]. (required)
   * - ``--genome``
     - string
     - (none)
     - Reference genome FASTA (run-all -g). (required)
   * - ``--te_library``
     - string
     - (none)
     - Transposon library FASTA (run-all -T). (required)
   * - ``--repeatmasker``
     - string
     - (none)
     - RepeatMasker .out for the reference (run-all --repeatmasker). Optional.
   * - ``--outdir``
     - string
     - (none)
     - Output directory; each sample is written to <outdir>/<sample>. (required)

Calling options
~~~~~~~~~~~~~~~

TE search, genome placement and insertion calling. Same names and defaults as relocaTE3 run-all.

Defaults reproduce RelocaTE2. Each option maps to the run-all flag named in its description.

.. list-table::
   :header-rows: 1
   :widths: 22 10 16 52

   * - Parameter
     - Type
     - Default
     - Description
   * - ``--te_aligner``
     - string
     - ``blat``
     - Aligner for the TE-library search (--te-aligner). One of: ``minimap2``, ``bwa``, ``bwamem2``, ``bwaaln``, ``bowtie2``, ``blat``.
   * - ``--genome_aligner``
     - string
     - ``bwaaln``
     - Aligner for genome re-alignment (--genome-aligner). blat is not supported here. One of: ``minimap2``, ``bwa``, ``bwamem2``, ``bwaaln``, ``bowtie2``.
   * - ``--te_opts``
     - string
     - ``''``
     - Extra options passed verbatim to the TE-search aligner (--te-opts), e.g. "-minIdentity=80" for blat.
   * - ``--genome_opts``
     - string
     - ``''``
     - Extra options passed verbatim to the genome aligner (--genome-opts), e.g. "-n 0.10" for bwa aln.
   * - ``--mismatch``
     - integer
     - ``2``
     - Allowed mismatches, TE alignment and read/genome comparison (--mismatch).
   * - ``--min_match``
     - integer
     - ``10``
     - Minimum alignment match length to the TE (--min-match).
   * - ``--min_trimmed``
     - integer
     - ``10``
     - Minimum trimmed flanking sequence length to retain (--min-trimmed).
   * - ``--size``
     - integer
     - ``500``
     - Library insert size for mate-pair-only insertions (--size).
   * - ``--min_mapq``
     - integer
     - ``0``
     - Minimum MAPQ for a read to be admitted as evidence (--min-mapq). 0 = no gate.
   * - ``--tsd``
     - string
     - ``UNK``
     - TSD mode (--tsd). Only UNK is supported: each TSD is inferred from the reads, as RelocaTE2 does. One of: ``UNK``, ``UKN``, ``unknown``.
   * - ``--require_both_junctions``
     - boolean
     - ``true``
     - Require junction reads on both sides of an insertion (--require-both-junctions).
   * - ``--target``
     - string
     - ``ALL``
     - Chromosome to analyze, or ALL (--target).
   * - ``--distance``
     - integer
     - ``3``
     - Reference-TE boundary window for dropping one-sided calls (--distance).
   * - ``--te_name``
     - string
     - ``repeat``
     - TE label used in output file names (--te-name).

BLAT scatter options
~~~~~~~~~~~~~~~~~~~~

How the BLAT TE search (te_aligner blat) is spread over tasks. Outputs are identical with any of these settings.

With blat_scatter, each sample's reads are split into chunk FASTAs (RELOCATE3_TE_SPLIT), the chunks are searched blat_chunks_per_task at a time (RELOCATE3_TE_SEARCH_CHUNK), and one task merges and trims them (RELOCATE3_TE_MERGE_TRIM), writing the files of the single RELOCATE3_TE_SEARCH task. Other TE aligners always run RELOCATE3_TE_SEARCH.

.. list-table::
   :header-rows: 1
   :widths: 22 10 16 52

   * - Parameter
     - Type
     - Default
     - Description
   * - ``--blat_scatter``
     - boolean
     - ``true``
     - Spread the BLAT TE search over many tasks. false: one RELOCATE3_TE_SEARCH task per sample, as relocaTE3 run.
   * - ``--blat_chunks_per_task``
     - integer
     - ``8``
     - Chunk FASTAs (of blat_chunk_reads reads) searched per RELOCATE3_TE_SEARCH_CHUNK task, 4 at once.
   * - ``--blat_chunk_reads``
     - integer
     - ``200000``
     - Reads per chunk FASTA (RelocaTE2's 200,000). Smaller chunks are for tests. (hidden from ``--help``)

Genotyping options
~~~~~~~~~~~~~~~~~~

Classify each insertion as homozygous, heterozygous or somatic.

On by default in the pipeline (opt-in with --genotype for run-all), so the cohort genotype matrix is filled. A sample's genotype_bam column reuses an existing alignment instead of realigning its reads, and such a sample is always genotyped, as run-all does with --genotype-bam.

.. list-table::
   :header-rows: 1
   :widths: 22 10 16 52

   * - Parameter
     - Type
     - Default
     - Description
   * - ``--genotype``
     - boolean
     - ``true``
     - Genotype each insertion (run-all --genotype). With false, only samples with a genotype_bam are genotyped.
   * - ``--excision``
     - boolean
     - ``false``
     - Also search for excision events that leave a footprint (--excision).
   * - ``--bcftools``
     - string
     - ``bcftools``
     - bcftools executable, used only with --excision (run-all --bcftools).

Report options
~~~~~~~~~~~~~~

The cohort tables and the MultiQC summary report written after all samples.

.. list-table::
   :header-rows: 1
   :widths: 22 10 16 52

   * - Parameter
     - Type
     - Default
     - Description
   * - ``--skip_multiqc``
     - boolean
     - ``false``
     - Do not write the MultiQC report (<outdir>/multiqc/). The cohort tables (<outdir>/cohort/) are still written.
   * - ``--multiqc_title``
     - string
     - (none)
     - Title of the MultiQC report. Default: RelocaTE3 (assets/multiqc_config.yml).

Container options
~~~~~~~~~~~~~~~~~

The image the docker, apptainer and singularity profiles run tasks in.

Tasks run in <container_registry>:<container_tag>. CI publishes dev (main branch), pr-<n>, sha-<short> and release versions (e.g. 0.3.0). Pin a reproducible run with a version or sha-<short> tag.

.. list-table::
   :header-rows: 1
   :widths: 22 10 16 52

   * - Parameter
     - Type
     - Default
     - Description
   * - ``--container_registry``
     - string
     - ``ghcr.io/stajichlab/relocate3``
     - Image repository, without the tag.
   * - ``--container_tag``
     - string
     - ``dev``
     - Image tag: dev, pr-<n>, sha-<short> or a release version.

Profile options
~~~~~~~~~~~~~~~

Settings used by specific profiles.

The pixi profile runs tasks without a container, with this environment's bin/ first on PATH.

.. list-table::
   :header-rows: 1
   :widths: 22 10 16 52

   * - Parameter
     - Type
     - Default
     - Description
   * - ``--pixi_env``
     - string
     - (none)
     - Pixi environment providing relocaTE3, the aligners and MultiQC (pixi profile). Default: <projectDir>/.pixi/envs/pipeline. (hidden from ``--help``)
   * - ``--test_data``
     - string
     - (none)
     - Test data directory (test profile: <projectDir>/tests/data). The test sample sheet refers to it as ${params.test_data}. (hidden from ``--help``)
   * - ``--slurm_partition``
     - string
     - ``intel``
     - SLURM partition tasks are submitted to (hpcc profile). conf/hpcc.config's resource limits fit intel nodes.

Generic options
~~~~~~~~~~~~~~~

Options that are not about the analysis.

.. list-table::
   :header-rows: 1
   :widths: 22 10 16 52

   * - Parameter
     - Type
     - Default
     - Description
   * - ``--help``
     - boolean or string
     - ``false``
     - Print the parameters and exit. With a parameter name (--help genome), print that parameter in detail. (hidden from ``--help``)

.. END GENERATED PARAMETERS


Profiles
--------

Choose exactly one software profile, and add ``hpcc`` or ``test`` as needed.

.. list-table::
   :header-rows: 1
   :widths: 15 85

   * - Profile
     - What it does
   * - ``apptainer``
     - Every task runs in ``<container_registry>:<container_tag>`` with
       Apptainer. Images are cached in ``NXF_APPTAINER_CACHEDIR`` (else
       ``<work-dir>/apptainer``). The pull may take up to 60 min.
   * - ``singularity``
     - As ``apptainer``, with Singularity. The cache variable is
       ``NXF_SINGULARITY_CACHEDIR``.
   * - ``docker``
     - As ``apptainer``, with Docker. Tasks run as your user, so output files
       are yours.
   * - ``pixi``
     - No container. Tasks run on the host with ``--pixi_env``'s ``bin/``
       first on ``PATH``. The default is the clone's ``.pixi/envs/pipeline``
       (``pixi install -e pipeline``). For development.
   * - ``hpcc``
     - UCR HPCC: each task is a SLURM job in ``--slurm_partition`` (default
       ``intel``). At most 50 jobs queued at once, submitted at most 10 per
       minute. Task resources are capped at 64 CPUs, 430 GB and 30 days.
       Combine with a container profile.
   * - ``test``
     - The Chr3 test data from ``tests/data`` and
       ``--outdir results/nextflow_test``. Caps tasks at 4 CPUs and 8 GB.
       Works from any launch directory.

**Profile order matters.** Profiles apply in command-line order, and for a
setting two profiles both set, the last one wins. ``test`` and ``hpcc`` both
set task resource caps:

- ``-profile test,hpcc,apptainer`` keeps hpcc's caps (64 CPUs, 430 GB);
- ``-profile hpcc,test,apptainer`` caps tasks at 4 CPUs and 8 GB.

Resources come from labels in ``conf/base.config``. A task killed for memory
or time (exit 137 or 140) is retried up to twice, with its memory and time
multiplied by the attempt number. Any other failure stops the run.


UCR HPCC launch
---------------

Run Nextflow itself in a small sbatch "head" job, never on the login node.
The head job submits one SLURM job per task and waits for them.

``run_relocate3_nf.sh``:

.. code-block:: bash

   #!/usr/bin/bash -l
   #SBATCH -p intel
   #SBATCH --cpus-per-task=2
   #SBATCH --mem=8G
   #SBATCH --time=1-00:00:00
   #SBATCH -o logs/relocate3_nf.%j.log
   set -euo pipefail
   date

   module load nextflow/26.04.6 apptainer/1.5.3

   # One shared image cache on bigdata: home quotas are small, and every
   # node must see the image.
   export NXF_APPTAINER_CACHEDIR=/bigdata/<lab>/<user>/apptainer-cache
   mkdir -p "${NXF_APPTAINER_CACHEDIR}" logs

   nextflow run stajichlab/RelocaTE3 -r v0.3.0rc1 -profile hpcc,apptainer \
     --container_tag 0.3.0rc1 \
     --slurm_partition intel \
     --input samples.csv \
     --genome reference.fa \
     --te_library mping.fa --te_name mping \
     --repeatmasker reference.fa.RepeatMasker.out \
     --outdir results \
     -work-dir /bigdata/<lab>/<user>/relocate3-work \
     -resume

Submit it from the directory that holds ``samples.csv``:

.. code-block:: bash

   mkdir -p logs && sbatch run_relocate3_nf.sh

Create ``logs/`` before submitting: SLURM opens the ``#SBATCH -o`` log before
the script runs, so with no ``logs/`` the job fails without writing a log.

Notes:

- ``NXF_APPTAINER_CACHEDIR``: Nextflow pulls the image once, on the head
  job's node, into this directory. Without it the image goes into
  ``<work-dir>/apptainer``, once per work directory.
- ``--slurm_partition`` sets the partition of the task jobs. The head job's
  partition is its own ``#SBATCH -p``. The resource caps of ``hpcc`` fit
  ``intel`` nodes.
- ``-work-dir`` holds every task's files. Put it on bigdata or scratch with
  room to spare (see `BLAT scatter`_ for the disk cost).
- ``-resume`` reuses finished tasks when you submit the same script again,
  for example after a failure or with more samples. It is harmless on the
  first run. Launch from the same directory: Nextflow keeps its run history
  in ``.nextflow/`` there.
- When the head job stops (for example at its time limit), Nextflow normally
  cancels the task jobs it submitted, when it is stopped cleanly (``scancel``
  of the head job sends SIGTERM). If it is killed abruptly, cancel leftover
  task jobs with ``scancel``. Submit the script again; ``-resume`` picks up
  from the finished tasks.
- The head job only submits and waits for the task jobs, so it needs little
  CPU and memory. Its 1-day limit fits the measured runs below; give a large
  cohort more with ``sbatch --time=3-00:00:00 run_relocate3_nf.sh``. Keep the
  limit no longer than needed: SLURM holds a job whose limit overlaps a
  scheduled maintenance reservation until the reservation ends.
- ``conf/hpcc.config``'s ``resourceLimits`` (64 CPUs, 430 GB, 30 days) is a
  ceiling, not a request: Nextflow lowers any larger request, including a
  retry's, to it. Each task asks for its label's resources from
  ``conf/base.config``.

Measured: 10 real rice samples with ``-profile hpcc,apptainer`` and image
``sha-2c9a448`` (head job 31830681, ``validation/real_rice/run_pipeline.sh``)
ran for 1 h 04 min (time spent pending in the queue is not counted). Every
sample's files were identical to the v0.2.0rc2 ``run-all`` results
(``validation/real_rice/compare_pipeline_to_rc2.sh``, job 31833986).


Outputs
-------

.. code-block:: text

   <outdir>/
     <sample>/            one directory per sample: the same files as run-all -o <sample>
     cohort/              tables over all genotyped samples
     multiqc/             relocate3_report.html and its data directory
     pipeline_info/       versions, environment and Nextflow's run reports

``<sample>/``
   The same layout and file names as ``relocaTE3 run-all`` (see
   :doc:`output` and the README's "Outputs" section). Nothing is added.
   Chunk files of the BLAT scatter, merged lanes and the genome index are not
   published.

``cohort/``
   The ``relocaTE3 cohort`` outputs over every characterized sample:
   ``cohort.characTErized.txt`` (all rows), ``cohort.genotype_matrix.tsv``
   (sites x samples) and ``cohort.summary_mqc.tsv`` (calls per sample and per
   genotype class). Written only when at least one sample was genotyped. See
   the README's "Cohort outputs" section.

``multiqc/relocate3_report.html``
   One summary report: calls per sample and per genotype class, per-sample
   counts (TE-containing reads, flanking reads, calls per output tier,
   reference calls) and software versions. The per-sample count table is also
   in ``multiqc/relocate3_report_data/multiqc_relocate3_sample_stats_table.txt``.
   ``--skip_multiqc`` turns it off; ``--multiqc_title`` sets its title.

``pipeline_info/``
   - ``software_versions.yml``: the version of every tool each process used,
     and of the pipeline and Nextflow.
   - ``software_environment.txt``: what ran the tasks: date, run name,
     pipeline version, revision and commit, profile, container engine and the
     image reference (or the pixi environment).
   - ``execution_report.html``, ``execution_timeline.html``,
     ``execution_trace.txt`` (per-task time, CPU and memory) and
     ``pipeline_dag.html``.

   These have fixed names and are overwritten by the next run into the same
   ``--outdir``. Copy them first if you need to keep them.


BLAT scatter
------------

With ``--te_aligner blat`` (the default), the TE search of one sample is
split over many tasks, so it no longer runs as one long job.

1. ``RELOCATE3_TE_SPLIT`` cuts each read side into chunks of 200,000 reads
   (``relocaTE3 te-split``). These are the chunks ``relocaTE3 run`` makes
   internally.
2. ``RELOCATE3_TE_SEARCH_CHUNK`` searches ``--blat_chunks_per_task`` chunks
   (default 8) per task, 4 at a time (``relocaTE3 te-search-chunk``).
3. ``RELOCATE3_TE_MERGE_TRIM`` merges each side's results in order and trims
   them, once per sample (``relocaTE3 te-merge-trim``).

The outputs are identical, file for file, to the single
``RELOCATE3_TE_SEARCH`` task. For one 30x rice sample with the riceTElib
library, the TE search took 13 h 13 min as one task and 1 h 21 min
scattered, with identical outputs (``docs/2026-11-05-nextflow-m4.md``).

- ``--blat_chunks_per_task N``: more chunks per task means fewer, longer
  jobs. It does not change the results. At 8, the 30x tasks took 33 to 42
  minutes each.
- **Disk.** The chunk files stay in the work directory until you clean it
  (``nextflow clean``, or delete ``-work-dir`` after the run). For one 30x
  rice sample (374 chunks) that is about 15 GB of chunk FASTAs plus 6 GB of
  chunk BAMs, about 21 GB in total. Leave room for this for every sample in
  flight.
- **When to turn it off** (``--blat_scatter false``, one
  ``RELOCATE3_TE_SEARCH`` task per sample): on a single machine without a
  scheduler, or on a cluster that limits the number of jobs, where many
  small jobs cost more than they save.
- Other TE aligners (``--te_aligner minimap2``, ``bwaaln``, ...) always run
  one ``RELOCATE3_TE_SEARCH`` task.


Reproducibility
---------------

Pin both the pipeline code and the image:

- ``-r v0.3.0rc1`` pins the pipeline (a tag or a commit). Without ``-r``,
  Nextflow uses whatever revision it pulled last.
- ``--container_tag 0.3.0rc1`` pins the image. A release tag's default
  ``--container_tag`` is its own version, but passing it makes the command
  self-describing.

Image tags of ``ghcr.io/stajichlab/relocate3``:

.. list-table::
   :header-rows: 1
   :widths: 20 80

   * - Tag
     - Meaning
   * - ``0.3.0rc1`` (a version)
     - Built from the release tag ``v0.3.0rc1``. Fixed.
   * - ``sha-<short>``
     - Built from one commit. Fixed.
   * - ``dev``
     - Latest ``main``. **Moves** with every push.
   * - ``pr-<n>``
     - Latest commit of pull request ``<n>``. **Moves** with every push.

Use a version or ``sha-<short>`` tag for results you keep. A moving tag can
also mix images within one analysis: Nextflow's ``-resume`` keys tasks on the
image *name*, so tasks cached from an older ``dev`` image are reused after
``dev`` moved.

What a run records, in ``<outdir>/pipeline_info/``:

- ``software_versions.yml``: tool versions per process;
- ``software_environment.txt``: revision, commit, profile and image
  reference. Nextflow does not record the image digest; the file says how to
  get it (``apptainer inspect`` of the cached image, or
  ``docker image inspect``).

The result tables are reproducible: the same inputs, version and settings
give the same tables. Intermediate BAMs are written by multi-threaded
aligners and can differ byte for byte between runs.


Troubleshooting
---------------

**Nextflow refuses to start: the version is too old.**
The pipeline needs Nextflow 26.04 or newer (``manifest.nextflowVersion``,
nf-schema 3.0.0). ``module load nextflow/26.04.6``, or use the clone's
``pixi run -e nf nextflow``.

**"Cannot cast object 'ghcr.io/stajichlab/relocate3:...' with class 'java.lang.String' to class 'java.util.Map'".**
A pipeline revision older than the nf-schema 3.0.0 upgrade (November 2026),
run with a container profile. Use a current revision (``-r v0.3.0rc1`` or
``-r main``; ``nextflow pull stajichlab/RelocaTE3`` refreshes a cached copy).

**"Missing required parameter(s)" or a schema error before any task.**
``--input``, ``--genome``, ``--te_library`` and ``--outdir`` are required, and
every file must exist. Values are checked against ``nextflow_schema.json``:
``--te_name`` and ``--target`` may contain only ``[A-Za-z0-9._-]``; ``--tsd``
takes only ``UNK`` (RelocaTE2 always infers each TSD from the reads).

**A sample-sheet file "does not exist", though it does.**
Relative paths in the sheet are resolved against the launch directory, not
the sheet's directory. Use absolute paths.

**"Sample '<name>': genotype_bam <file> has no index".**
Create the index beside the file: ``samtools index <file>``. Only
``<file>.bai``, ``<file>.csi`` and ``<file>.crai`` are found; ``<stem>.bai``
(without ``.bam``) is not.

**"Sample name must be provided, start with a letter or digit ...".**
Rename the sample: letters, digits, ``.``, ``_`` and ``-`` only, not
``cohort``, ``multiqc`` or ``pipeline_info``.

**"Sample '<name>': its sample-sheet rows mix paired-end and single-end reads",
or "... give different genotype_bam values".**
All rows (lanes) of one sample must agree. Fix the sheet.

**"-profile pixi: no relocaTE3 at ...".**
Run ``pixi install -e pipeline`` in the clone, or point ``--pixi_env`` at an
environment that has ``relocaTE3``. The ``pixi`` profile works only from a
clone, not with ``nextflow run stajichlab/RelocaTE3``.

**``--genotype false`` (or another ``false``) seems to be ignored.**
Pipeline revisions before the Nextflow M4 fix (November 2026) read
``--<flag> false`` as true on Nextflow 26.04+. Current revisions type these
settings, so ``false`` works and a word other than true/false stops the run.
Note that a sample with a ``genotype_bam`` is always genotyped.

**The image pull is slow, times out, or fills your home directory.**
Set ``NXF_APPTAINER_CACHEDIR`` (or ``NXF_SINGULARITY_CACHEDIR``) to a shared
directory with space. The profiles allow 60 minutes per pull. On UCR HPCC,
``mkdir -p logs && sbatch containers/pull_and_check.sh <tag>`` (from a clone)
pulls the image
into that cache ahead of a run and checks the tools in it.

**A task was killed for memory or time.**
Exit 137 or 140 is retried twice with more memory and time. If it still
fails, raise the label's resources in a small config file and pass it with
``-c``:

.. code-block:: groovy

   // more.config
   process {
       withName: 'RELOCATE3_TE_SEARCH' {
           memory = { 32.GB * task.attempt }
           time = { 48.h * task.attempt }
       }
   }

.. code-block:: bash

   nextflow -c more.config run stajichlab/RelocaTE3 -r v0.3.0rc1 \
     -profile hpcc,apptainer --container_tag 0.3.0rc1 \
     --input samples.csv --genome reference.fa --te_library mping.fa \
     --outdir results -resume

**``RELOCATE3_TE_SPLIT`` fails for one sample.**
A sample with no reads fails here instead of silently giving no output.
``run-all`` fails on such a sample too. Check its FASTQ files.

**The work directory fills the disk.**
The BLAT scatter keeps about 21 GB per 30x sample in the work directory.
Run fewer samples at once, use ``--blat_scatter false``, or clean finished
runs with ``nextflow clean``.

**``--help`` starts with two stray description lines.**
Cosmetic, from nf-schema 3.0.0's help printer. The parameter list below them
is complete.

**The previous run's report is gone.**
``pipeline_info/`` files are overwritten by the next run into the same
``--outdir``.


For agents and scripts
----------------------

- **Parameters.** ``nextflow_schema.json`` (JSON Schema, draft 2020-12) is the
  full list, with types, defaults, allowed values and patterns.
  ``nextflow run stajichlab/RelocaTE3 -r <rev> --help`` prints it, and
  ``--help <parameter>`` prints one. The sample-sheet schema is
  ``assets/schema_input.json``.
- **Params file.** Put the settings in a YAML (or JSON) file and pass it with
  ``-params-file``. Values there are typed, so booleans are ``true`` /
  ``false``:

  .. code-block:: yaml

     input: /data/project/samples.csv
     genome: /data/ref/reference.fa
     te_library: /data/ref/mping.fa
     te_name: mping
     repeatmasker: /data/ref/reference.fa.RepeatMasker.out
     outdir: /data/project/results
     genotype: true
     container_tag: 0.3.0rc1

  .. code-block:: bash

     nextflow run stajichlab/RelocaTE3 -r v0.3.0rc1 -profile hpcc,apptainer \
       -params-file params.yaml -resume

- **Check before running.** ``-preview`` validates the parameters and the
  sample sheet and exits without running any task.
- **Provenance.** ``<outdir>/pipeline_info/software_environment.txt`` (revision,
  commit, profile, image), ``software_versions.yml`` (tool versions),
  ``execution_trace.txt`` (one line per task: status, exit code, time,
  memory). ``nextflow log`` in the launch directory lists past runs.
- **Exit status.** ``nextflow run`` exits non-zero when any task fails for
  good. A run that finished has every sample's
  ``<outdir>/<sample>/results/`` and, when samples were genotyped,
  ``<outdir>/cohort/``.
