Python API Reference
====================

.. toctree::
   :maxdepth: 1
   :caption: Modules:
   :hidden:

   uffpsim_database
   uffpsim_search_engine
   uffpsim_fingerprints
   uffpsim_fp_supplier
   uffpsim_fp_parallel_supplier
   uffpsim_web_app

Summary Table
-------------

.. list-table::
   :widths: 25 75
   :header-rows: 1

   * - Module
     - Description
   * - :mod:`uffpsim` (package)
     - Package-level exports for quick imports.
   * - :doc:`Database <uffpsim_database>`
     - Functions for creating, updating, and querying HDF5 databases of molecular fingerprints.
   * - :doc:`Search Engine <uffpsim_search_engine>`
     - Class for performing single and batch similarity searches on uffpsim databases.
   * - :doc:`Fingerprints <uffpsim_fingerprints>`
     - Fingerprint calculator supporting Morgan, RDKit, AtomPair, Torsion, MACCSKeys, Avalon, and Pattern fingerprint types.
   * - :doc:`FP Supplier <uffpsim_fp_supplier>`
     - Serial fingerprint generation from SDF, SMILES, InChI files or Python iterables.
   * - :doc:`Parallel FP Supplier <uffpsim_fp_parallel_supplier>`
     - Parallel fingerprint generation using multiple Python processes via ``ProcessPoolExecutor``.
   * - :doc:`Web App <uffpsim_web_app>`
     - Flask-based web application for interactive search.

Quick Start Example
-------------------

.. code-block:: python

    from uffpsim import (
        create_database,
        create_database_parallel,
        redo_inner_clustering,
        get_database_info,
        UFFPSimSearchEngine,
        FPCalculator,
        FP_TYPE_DEFAULT_PARAMETERS,
    )
    from uffpsim.fingerprints import FP_TYPES

    # 1. Create database
    create_database(
        input_file='molecules.sdf.gz',
        db_file='mydb.h5',
        fp_type='Morgan',
        fp_params={'fpSize': 2048},
    )

    # 2. Query database
    engine = UFFPSimSearchEngine('mydb.h5')
    hits = engine.search("CC(=O)Oc1ccccc1C(=O)O", threshold=0.7, limit_by=5)
    for compound_id, score in hits:
        print(f"{compound_id}: {score:.4f}")
