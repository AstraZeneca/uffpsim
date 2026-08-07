Introduction
============

Searching for chemically similar compounds is a fundamental task in computer-aided drug discovery, 
used throughout hit identification, hit expansion and lead optimisation.
With the growth of public databases such as ChEMBL, PubChem and ZINC into hundreds of millions or billions of molecules — 
and resources such as Enamine REAL exceeding 13 billion compounds — exact fingerprint similarity search remains both indispensable and increasingly expensive.
The BitBound algorithm (Swamidass and Baldi, 2007) addresses this by pruning candidates using bounds on the maximum achievable Tanimoto similarity, and underpins widely used libraries such as FPSim2 and chemfp.
Hardware acceleration (SIMD, GPUs) further speeds up individual comparisons, but does not reduce the amount of work that grows with database size.

**UFFPSim** (Ultrafast Fingerprint Similarity) extends BitBound with a second pruning stage based on clustered fingerprints within each popcount bin, discarding large regions of the search space before exact Tanimoto evaluation.
It offers both in-memory search for high-throughput screening and disk-based search for databases that exceed available RAM, enabling exact queries on databases with up to a billion compounds via an HDF5 store, a Python API and a CLI.
Benchmarks against FPSim2 on ChEMBL, PubChem and a one-billion-compound ZINC20 subset show consistently lower query latency while preserving exact results; optional AVX512-VPOPCNTDQ support is available but complementary to the core algorithm.

Architecture Overview
---------------------

uffpsim consists of three main components:

1. **Fingerprint Calculation** (``uffpsim.fingerprints``, ``uffpsim.fp_supplier``, ``uffpsim.fp_parallel_supplier``)
   - Computes molecular fingerprints using RDKit and Avalon libraries
   - Supports serial and parallel fingerprint generation from SDF, SMILES, or InChI files

2. **Database Management** (``uffpsim.database``)
   - Creates and manages HDF5 databases of molecular fingerprints
   - Performs inner clustering to group similar molecules by precomputed pop-count bins
   - Supports incremental database updates

3. **Search Engine** (``uffpsim.search_engine``)
   - Loads databases and performs single or batch similarity searches
   - Returns top-k matches above a given Tanimoto similarity threshold
   - Available as both Python API and CLI tool

Supported Fingerprint Types
----------------------------

The following fingerprint types are supported:

+------------------------+---------------------------------------------------+
| Type                   | Description                                       |
+========================+===================================================+
| Morgan                 | Circular fingerprints (ECFP-style) [default]      |
+------------------------+---------------------------------------------------+
| RDKit                  | Path-based fingerprints                           |
+------------------------+---------------------------------------------------+
| AtomPair               | Pairwise atom descriptor fingerprints             |
+------------------------+---------------------------------------------------+
| TopologicalTorsion     | Extended connectivity torsion fingerprints        |
+------------------------+---------------------------------------------------+
| MACCSKeys              | 166-bit MACCS structural keys                     |
+------------------------+---------------------------------------------------+
| Avalon                 | Avalon topological fingerprints                   |
+------------------------+---------------------------------------------------+
| RDKPatternFingerprint  | Pattern-based fingerprints                        |
+------------------------+---------------------------------------------------+

Input Formats
-------------

uffpsim supports the following molecular input formats:

* **SDF / SDF.GZ**: Structured Data File format (with optional gzip compression)
* **SMILES / SMILES with IDs**: One molecule per line, optionally with compound IDs as the second column
* **InChI**: One InChI string per line (optionally with compound IDs as the second column)

Database Options
----------------

When creating a database, users can control:

* **Fingerprint type and parameters** (e.g., Morgan radius, fingerprint size)
* **Molecule ID handling** (auto-generated serial numbers or from a property field)
* **Max molecule ID length** (`mol_id_max_chars`) — must be a multiple of 8 minus one (7, 15, 23, 31...)
* **Inner clustering threshold** — molecules with Tanimoto similarity above this threshold are grouped together for faster search
* **Clustering mode** — `memory` stores clustering data in RAM; `disk` persists it to the HDF5 file
* **Parallel processing** — distributed fingerprint generation via Python ``ProcessPoolExecutor`` and OpenMP-based parallel clustering

Command-Line Interface
----------------------

+--------------------------------------+--------------------------------------------------------------------+
| Command                              | Description                                                        |
+======================================+====================================================================+
| :doc:`cli/create-database`           | Create a new uffpsim HDF5 database from a molecular input file.    |
+--------------------------------------+--------------------------------------------------------------------+
| :doc:`cli/redo-clustering`           | Redo inner clustering for an existing database with a new          |
|                                      | similarity threshold.                                              |
+--------------------------------------+--------------------------------------------------------------------+
| :doc:`cli/build-mol-id-index-table`  | Build serialized MolIdIndexTable for an existing database.         |
+--------------------------------------+--------------------------------------------------------------------+
| :doc:`cli/search`                    | Search one or more SMILES against the database and write hits      |
|                                      | to a CSV file.                                                     |
+--------------------------------------+--------------------------------------------------------------------+
| :doc:`cli/launch-web-app`            | Launch the uffpsim interactive web application for visual          |
|                                      | similarity searching.                                              |
+--------------------------------------+--------------------------------------------------------------------+


Web Application
---------------

A built-in Flask web application provides an interactive interface for similarity searching. See :doc:`demo_webapp` for a guide to running it locally.

Python API Quick Reference
--------------------------

.. code-block:: python

    from uffpsim import create_database, UFFPSimSearchEngine

    # Create database
    create_database(
        input_file='chembl_33.sdf.gz',
        db_file='chembl.h5',
        fp_type='Morgan',
        fp_params={'fpSize': 2048, 'radius': 2},
    )

    # Search
    engine = UFFPSimSearchEngine('chembl.h5')
    results = engine.search("Cc1cc(-n2ncc(=O)[nH]c2=O)ccc1C(=O)c1ccc(C#N)cc1", 0.8, limit_by=10)
