Examples
========

This page contains practical examples for creating databases and performing searches with uffpsim.

Downloading Example Data
------------------------

The following command downloads the ChEMBL 33 database to use as an example:

.. code-block:: bash

    curl https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases/chembl_33/chembl_33.sdf.gz -o ./chembl_33.sdf.gz
    ls -artl chembl_33.sdf.gz
    # -rw-rw-r-- 1 user user 770312723 Jun 11 12:28 chembl_33.sdf.gz

Molecule ID Configuration
-------------------------

IDs are created automatically as string serial numbers starting from ``'1'``. The parameter ``mol_id_max_chars`` sets the maximum number of characters in a molecule ID and must be a multiple of 8 minus one (e.g., 7, 15, 23, 31). This value is fixed at database creation time and cannot be changed afterwards.

Molecule IDs are stored inline in the same fixed-size array as the fingerprint bits, so every entry in the database reserves exactly ``mol_id_max_chars + 1`` bytes for the ID regardless of actual ID length. Choose a value large enough to fit your longest molecule ID -- for example, ChEMBL IDs (``CHEMBL1234567``, 13 chars) fit within ``mol_id_max_chars=15``.

.. warning::
    IDs longer than ``mol_id_max_chars`` are silently truncated, which can cause incorrect lookups. Choose a value large enough for your dataset upfront.

**Memory trade-off**: each additional 8 characters of headroom costs 8 bytes per molecule. As a ballpark example, for a 2M-molecule database, increasing from ``mol_id_max_chars=15`` to ``mol_id_max_chars=31`` adds ~16 MB to the database size.

Creating a Database
-------------------

Serial Fingerprints, Single-Threaded Clustering
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Fingerprints are calculated one molecule at a time and clustering runs in a single thread. This is the simplest option and uses the least memory.

.. code-block:: python

    from uffpsim import create_database

    create_database(
        input_file='chembl_33.sdf.gz',
        db_file='chembl_2048b.h5',
        fp_type='Morgan',
        fp_params={'fpSize': 2048, 'radius': 2},
        gen_ids=True,
        mol_id_prop=None,
        mol_id_max_chars=15,
        inner_clustering_threshold=0.15,
    )

Serial Fingerprints, Parallel Clustering
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Fingerprint calculation is still serial, but the inner-clustering step uses multiple OpenMP threads within the same process. Set ``OMP_NUM_THREADS`` to control the thread count, e.g., ``export OMP_NUM_THREADS=4``.

.. code-block:: python

    from uffpsim import create_database

    create_database(
        input_file='chembl_33.sdf.gz',
        db_file='chembl_2048b.h5',
        fp_type='Morgan',
        fp_params={'fpSize': 2048, 'radius': 2},
        gen_ids=True,
        mol_id_prop=None,
        mol_id_max_chars=15,
        inner_clustering_threshold=0.15,
        cluster_parallel=True,
    )

Parallel Fingerprints, Parallel Clustering
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Fingerprint calculation is distributed across ``workers`` independent Python processes using ``ProcessPoolExecutor``, which bypasses the GIL and speeds up the RDKit computation step. Writing to the HDF5 file and inner-clustering remain serial and parallel (via OpenMP) respectively. Note that spawning multiple processes requires more memory than the serial approach.

.. code-block:: python

    from uffpsim import create_database_parallel

    create_database_parallel(
        input_file='chembl_33.sdf.gz',
        db_file='chembl_2048b.h5',
        workers=4,
        fp_type='Morgan',
        fp_params={'fpSize': 2048, 'radius': 2},
        gen_ids=True,
        mol_id_prop=None,
        mol_id_max_chars=15,
        inner_clustering_threshold=0.15,
        cluster_parallel=True,
    )

Redoing Inner Clustering
--------------------------

Inner-clustering can be re-performed with a different threshold. It is useful to benchmark different threshold values for a new database. The clustering data is written back to the same file.

.. code-block:: python

    from uffpsim import redo_inner_clustering

    redo_inner_clustering("chembl_2048b.h5", 0.1, cluster_parallel=True)

Searching a Database
--------------------

Single-Query Search
~~~~~~~~~~~~~~~~~~~

Each time only one SMILES string can be used as input.

.. code-block:: python

    from uffpsim import UFFPSimSearchEngine

    search_engine = UFFPSimSearchEngine("chembl_2048b.h5")

    # Print few parameters from database
    print(search_engine.fp_store.fp_params_json)
    print(search_engine.fp_store.mol_id_max_chars)
    print(search_engine.fp_store.fp_bits_size)
    print(search_engine.fp_store.inner_clustering_threshold)

    # Similarity threshold of 0.8 and return only up to one hit
    result_1 = search_engine.search("Cc1cc(-n2ncc(=O)[nH]c2=O)ccc1C(=O)c1ccc(C#N)cc1", 0.8, limit_by=1)
    print(result_1)

    # Similarity threshold of 0.6 and return only up to 10 hits
    result_2 = search_engine.search("Cc1cc(-n2ncc(=O)[nH]c2=O)ccc1C(=O)c1ccc(C#N)cc1", 0.6, limit_by=10)
    print(result_2)

.. code-block:: python

    # Example of result from search
    >>> search_engine.search("Cc1cc(-n2ncc(=O)[nH]c2=O)ccc1C(=O)c1ccc(C#N)cc1", 0.8, limit_by=1)
    [('2', 1.0)]  # '2' matches the second SMILES string in `tests/data/10mols.smi`

Batch Search
~~~~~~~~~~~~

Searching in batch mode could potentially speed-up the search by up to two times depending on the nature of queries.

.. code-block:: python

    from uffpsim import UFFPSimSearchEngine

    search_engine = UFFPSimSearchEngine("chembl_2048b.h5")

    # List of SMILES to be searched
    smiles_list = [
        "Cc1cc(-n2ncc(=O)[nH]c2=O)ccc1C(=O)c1ccccc1Cl",
        "Cc1cc(-n2ncc(=O)[nH]c2=O)ccc1C(=O)c1ccc(C#N)cc1",
        "Cc1cc(-n2ncc(=O)[nH]c2=O)cc(C)c1C(O)c1ccc(Cl)cc1",
        "Cc1ccc(C(=O)c2ccc(-n3ncc(=O)[nH]c3=O)cc2)cc1",
        "Cc1cc(-n2ncc(=O)[nH]c2=O)ccc1C(=O)c1ccc(Cl)cc1",
        "Cc1cc(-n2ncc(=O)[nH]c2=O)ccc1C(=O)c1ccccc1",
        "Cc1cc(Br)ccc1C(=O)c1ccc(-n2ncc(=O)[nH]c2=O)cc1Cl",
        "O=C(c1ccc(Cl)cc1Cl)c1ccc(-n2ncc(=O)[nH]c2=O)cc1Cl",
        "CS(=O)(=O)c1ccc(C(=O)c2ccc(-n3ncc(=O)[nH]c3=O)cc2Cl)cc1",
        "c1cc2cc(c1)-c1cccc(c1)C[n+]1ccc(c3ccccc31)NCCCCCCCCCCNc1cc[n+](c3ccccc13)C2",
    ]

    # Similarity threshold of 0.8 and return only up to one hit
    result_1 = search_engine.batch_search(smiles_list, 0.8, limit_by=1)
    print(result_1)

    # Similarity threshold of 0.6 and return up to 10 hits
    result_2 = search_engine.batch_search(smiles_list, 0.6, limit_by=10)
    print(result_2)
