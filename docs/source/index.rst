RelocaTE3
=========

RelocaTE3 identifies transposable element (TE) insertion polymorphisms from
short-read resequencing data at single-base resolution.  It is a modern
pure-Python reimplementation of `RelocaTE2
<https://github.com/JinfengChen/RelocaTE2>`_ that requires only **minimap2**
and **samtools** as external tools — no blat, bwa, bowtie2, seqtk, or Perl.

**Status:** reference implementation.  On the rice Chr3 2 Mb benchmark the
tool recovers ~178/200 simulated insertions (~89 % recall, ~90 % precision);
see ``tests/acceptance_test.py``.

.. toctree::
   :maxdepth: 2
   :caption: Contents

   installation
   usage
   output
   migration
   api

Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
