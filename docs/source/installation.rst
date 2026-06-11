Installation
============

Prerequisites
-------------

RelocaTE3 requires **Python ≥ 3.10** and two external bioinformatics tools:

- `minimap2 <https://github.com/lh3/minimap2>`_ ≥ 2.24
- `samtools <https://www.htslib.org/>`_ ≥ 1.17

The optional false-positive filtering and BED output also uses:

- `bedtools <https://bedtools.readthedocs.io/>`_ ≥ 2.30

With pixi (recommended)
-----------------------

`pixi <https://pixi.sh>`_ pins all tool versions automatically:

.. code-block:: bash

   git clone https://github.com/stajichlab/RelocaTE3
   cd RelocaTE3
   pixi install
   pixi run relocaTE3 --help

With pip
--------

If minimap2 and samtools are already on your ``PATH``:

.. code-block:: bash

   pip install -e .
   relocaTE3 --help

On an HPC cluster (module-based)
---------------------------------

.. code-block:: bash

   module load minimap2 samtools bedtools
   pip install --user -e .
