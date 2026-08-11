Usage
========

This page contains practical usage instructions for creating databases and performing searches with uffpsim.

.. contents:: Table of Contents
   :local:
   :depth: 2

Downloading Example Data
------------------------

The following command downloads the ChEMBL 35 database to use as an example:

.. tab:: Bash

   .. code-block:: bash

      curl https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases/chembl_35/chembl_35.sdf.gz -o ./chembl_35.sdf.gz
      ls -artl chembl_35.sdf.gz
      # -rw-rw-r-- 1 user user 770312723 Jun 11 12:28 chembl_35.sdf.gz

Molecule ID Configuration
-------------------------

IDs are created automatically as string serial numbers starting from ``'1'``. The parameter ``mol_id_max_chars`` 
sets the maximum number of characters in a molecule ID and must be a multiple of 8 minus one (e.g., 7, 15, 23, 31). 
This value is fixed at database creation time and cannot be changed afterwards.

Molecule IDs are stored inline in the same fixed-size array as the fingerprint bits, so every entry in the database 
reserves exactly ``mol_id_max_chars + 1`` bytes for the ID regardless of actual ID length. Choose a value large 
enough to fit your longest molecule ID -- for example, ChEMBL IDs (``CHEMBL1234567``, 13 chars)
fit within ``mol_id_max_chars=15``.

.. warning::
    IDs longer than ``mol_id_max_chars`` are silently truncated, which can cause incorrect lookups. 
    Choose a value large enough for your dataset upfront.

**Memory trade-off**: each additional 8 characters of headroom costs 8 bytes per molecule. 
As a ballpark example, for a 2M-molecule database, increasing from ``mol_id_max_chars=15`` 
to ``mol_id_max_chars=31`` adds ~16 MB to the database size.


Building Dataset 
------------------

Serial Fingerprints, Single-Threaded Clustering
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Fingerprints are calculated one molecule at a time and clustering runs in a single thread. This is the simplest option and uses the least memory.

.. tab:: Python

   .. code-block:: python

      from uffpsim import create_database

      create_database(
          input_file='chembl_35.sdf.gz',
          db_file='chembl_2048b.h5',
          fp_type='Morgan',
          fp_params={'fpSize': 2048, 'radius': 2},
          gen_ids=True,
          mol_id_prop=None,
          mol_id_max_chars=15,
          inner_clustering_threshold=0.15,
      )

.. tab:: Bash

   .. code-block:: bash

      uffpsim create-database \
          --input-file chembl_35.sdf.gz \
          --db-file chembl_2048b.h5 \
          --fp-type Morgan \
          --fp-size 2048 \
          --radius 2 \
          --mol-id-max-chars 15 \
          --cluster-threshold 0.15

Supported fingerprint types: ``Morgan``, ``RDKit``, ``AtomPair``, ``TopologicalTorsion``, ``MACCSKeys``, ``Avalon``, and ``RDKPatternFingerprint`` (default: ``Morgan``).

Serial Fingerprints, Parallel Clustering
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Fingerprint calculation is still serial, but the inner-clustering step uses multiple OpenMP threads within the same process. Set ``OMP_NUM_THREADS`` to control the thread count.

.. tab:: Python

   .. code-block:: python

      from uffpsim import create_database

      create_database(
          input_file='chembl_35.sdf.gz',
          db_file='chembl_2048b.h5',
          fp_type='Morgan',
          fp_params={'fpSize': 2048, 'radius': 2},
          gen_ids=True,
          mol_id_prop=None,
          mol_id_max_chars=15,
          inner_clustering_threshold=0.15,
          cluster_parallel=True,
      )

.. tab:: Bash

   .. code-block:: bash

      export OMP_NUM_THREADS=4
      uffpsim create-database \
          --input-file chembl_35.sdf.gz \
          --db-file chembl_2048b.h5 \
          --fp-type Morgan \
          --fp-size 2048 \
          --radius 2 \
          --mol-id-max-chars 15 \
          --cluster-threshold 0.15 \
          --cluster-parallel

Parallel Fingerprints, Parallel Clustering
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Fingerprint calculation is distributed across ``workers`` independent Python processes using ``ProcessPoolExecutor``, which bypasses the GIL and speeds up the RDKit computation step. Writing to the HDF5 file and inner-clustering remain serial and parallel (via OpenMP) respectively. Note that spawning multiple processes requires more memory than the serial approach.

.. tab:: Python

   .. code-block:: python

      from uffpsim import create_database_parallel

      create_database_parallel(
          input_file='chembl_35.sdf.gz',
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

.. tab:: Bash

   .. code-block:: bash

      export OMP_NUM_THREADS=4
      uffpsim create-database \
          --input-file chembl_35.sdf.gz \
          --db-file chembl_2048b.h5 \
          --fp-type Morgan \
          --fp-size 2048 \
          --radius 2 \
          --mol-id-max-chars 15 \
          --cluster-threshold 0.15 \
          --cluster-parallel \
          --workers 4

Parallel Fingerprints, Parallel Clustering (**CUDA**)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Inner clustering can be performed using the CUDA also. The CUDA implementation is faster than the OpenMP implementation.
The CUDA implementation is available only when UFFPSim is compiled with CUDA support.

.. tab:: Python

   .. code-block:: python

      from uffpsim import create_database_parallel

      create_database_parallel(
          input_file='chembl_35.sdf.gz',
          db_file='chembl_2048b.h5',
          workers=4,
          fp_type='Morgan',
          fp_params={'fpSize': 2048, 'radius': 2},
          gen_ids=True,
          mol_id_prop=None,
          mol_id_max_chars=15,
          inner_clustering_threshold=0.15,
          cluster_mode='cuda',
          cluster_parallel=True,
      )

.. tab:: Bash

   .. code-block:: bash

      export OMP_NUM_THREADS=4
      uffpsim create-database \
          --input-file chembl_35.sdf.gz \
          --db-file chembl_2048b.h5 \
          --fp-type Morgan \
          --fp-size 2048 \
          --radius 2 \
          --mol-id-max-chars 15 \
          --cluster-threshold 0.15 \
          --cluster-parallel \
          --cluster-mode cuda \
          --workers 4

Redoing Inner Clustering
------------------------

Inner-clustering can be re-performed with a different threshold. It is useful to benchmark different threshold 
values for a new database. The clustering data is written back to the same file.

.. note:: Three types of ``--cluster-mode`` are supported: ``memory``, ``disk``, and ``cuda``.

.. tab:: Python

   .. code-block:: python

      from uffpsim import redo_inner_clustering

      redo_inner_clustering("chembl_2048b.h5", 0.1, cluster_parallel=True)

.. tab:: Bash

   .. code-block:: bash

      uffpsim redo-clustering --db-file chembl_2048b.h5 --threshold 0.1 --cluster-parallel


Building a Molecule ID Index Table
----------------------------------

An index table maps molecule IDs to their row positions in the database, enabling faster lookups for SMILES.
It is required to enable fetching of SMILES of hit compounds, and therefore recommended if databse is to be used 
with the **web-app and REST API**. Build this table via CLI before using it from either interface:

.. tab:: Bash

   .. code-block:: bash

      uffpsim build-mol-id-index-table --db-file chembl_2048b.h5


.. warning:: This steps load the entire database into memory. For large databases, ensure sufficient RAM is available.
             Also, the index table persists in the database file, therefore it is not necessary to rebuild it for the 
             same database unless the database is modified. It also increases the database file size.

Searching a Database
--------------------

Single-Query Search
~~~~~~~~~~~~~~~~~~~

Each time only one SMILES string can be used as input. Use ``-o`` to write CLI results to a CSV file.
``mode`` option control whether the database is loaded into memory or read from disk on demand. 
The default is ``memory`` if supported by the database, otherwise ``disk``.

.. tab:: Python

   .. code-block:: python

      from uffpsim import UFFPSimSearchEngine

      search_engine = UFFPSimSearchEngine("chembl_2048b.h5", mode='memory')

      print(search_engine.fp_store.fp_params_json)
      print(search_engine.fp_store.mol_id_max_chars)
      print(search_engine.fp_store.fp_bits_size)
      print(search_engine.fp_store.inner_clustering_threshold)

      result_1 = search_engine.search("Cc1cc(-n2ncc(=O)[nH]c2=O)ccc1C(=O)c1ccc(C#N)cc1", 0.8, limit_by=1)
      print(result_1)

      result_2 = search_engine.search("Cc1cc(-n2ncc(=O)[nH]c2=O)ccc1C(=O)c1ccc(C#N)cc1", 0.6, limit_by=10)
      print(result_2)

.. tab:: Bash

   .. code-block:: bash

      uffpsim search \
          --db-file chembl_2048b.h5 \
          --query "Cc1cc(-n2ncc(=O)[nH]c2=O)ccc1C(=O)c1ccc(C#N)cc1" \
          --threshold 0.6 \
          --limit 10 \
          --output-csv results.csv

CSV output writes columns ``mol_id,similarity`` with one row per hit (no header line).


Batch Search
~~~~~~~~~~~~

Searching in batch mode could potentially speed-up the search by up to two times depending on the nature of queries. 
Provide a file with one SMILES per line using ``--input-file`` (CLI) or pass a list (API). 
CLI results are written to CSV via ``--output-csv``.

.. tab:: Python

   .. code-block:: python

      from uffpsim import UFFPSimSearchEngine

      search_engine = UFFPSimSearchEngine("chembl_2048b.h5", mode='memory')

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

      result_1 = search_engine.batch_search(smiles_list, 0.8, limit_by=1)
      print(result_1)

      result_2 = search_engine.batch_search(smiles_list, 0.6, limit_by=10)
      print(result_2)

.. tab:: Bash

   .. code-block:: bash

      cat > queries.smi << EOF
      Cc1cc(-n2ncc(=O)[nH]c2=O)ccc1C(=O)c1ccccc1Cl
      Cc1cc(-n2ncc(=O)[nH]c2=O)ccc1C(=O)c1ccc(C#N)cc1
      Cc1cc(-n2ncc(=O)[nH]c2=O)cc(C)c1C(O)c1ccc(Cl)cc1
      EOF

      uffpsim search \
          --db-file chembl_2048b.h5 \
          --input-file queries.smi \
          --threshold 0.6 \
          --limit 10 \
          --output-csv results.csv

Batch Search Results with SMILES of Hits
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. tab:: Python

   .. code-block:: python

      from uffpsim import UFFPSimSearchEngine

      search_engine = UFFPSimSearchEngine("chembl_2048b.h5", mode='memory')

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

      idToSmiles = {}
      result_1 = search_engine.batch_search(smiles_list, 0.8, limit_by=1)
      print(result_1)

      # Fetch SMILES for all hits in the results
      for r in result_1:
          for idn, _ in r:
              if idn not in idToSmiles:
                  idToSmiles[idn] = search_engine.get_smiles_for_id(idn)

      result_2 = search_engine.batch_search(smiles_list, 0.6, limit_by=10)
      print(result_2)

      # Fetch SMILES for all hits in the results
      for r in result_2:
          for idn, _ in r:
              if idn not in idToSmiles:
                  idToSmiles[idn] = search_engine.get_smiles_for_id(idn)

      print(idToSmiles)

.. tab:: Bash

   .. code-block:: bash

      cat > queries.smi << EOF
      Cc1cc(-n2ncc(=O)[nH]c2=O)ccc1C(=O)c1ccccc1Cl
      Cc1cc(-n2ncc(=O)[nH]c2=O)ccc1C(=O)c1ccc(C#N)cc1
      Cc1cc(-n2ncc(=O)[nH]c2=O)cc(C)c1C(O)c1ccc(Cl)cc1
      EOF

      uffpsim search \
          --db-file chembl_2048b.h5 \
          --input-file queries.smi \
          --threshold 0.6 \
          --limit 10 \
          --output-csv results.csv \
          --include-hit-smiles

Web-App and REST API
----------------------
A minimal web application is included for interactive exploration of a database.
It can be launched via CLI. The REST API endpoint is also available for programmatic access.
It enables users to perform similarity search and retrieve hit details including SMILES and images.

.. note:: It enables ultra-fast search because the database is loaded into memory once at the start of the web-app.

The web-app can be launched as follows:

.. tab:: Bash

   .. code-block:: bash

      uffpsim launch-web-app -d chembl_2048b.h5 --results full -p 5000


Open ``http://localhost:5000`` in your browser. Type a SMILES string, adjust the similarity slider, and click on hits to inspect structures.

Using REST API interface: ``http://localhost:5000/api/search`` endpoint can be used for programmatic access to the same search functionality.
Payload example is as follows:

.. code-block:: json

    {
        "threshold":0.6,
        "limit_by":10,
        "smiles_text":"CC(=O)Oc1ccccc1C(=O)O"
    }


Once the app is started, it also run a REST API server. The query can be sent via ``POST`` request to the endpoint using curl as follows:

.. tab:: Bash

   .. code-block:: bash

      curl -X POST http://localhost:5000/api/search \
          -H "Content-Type: application/json" \
          -d '{"threshold":0.6,"limit_by":10,"smiles_text":"CC(=O)Oc1ccccc1C(=O)O"}'

For multiple queries, the payload can be modified as follows by adding newline tag separated SMILES strings in the ``smiles_text`` field:

.. tab:: Bash

   .. code-block:: bash

      curl -X POST http://localhost:5000/api/search \
          -H "Content-Type: application/json" \
          -d '{"threshold":0.6,"limit_by":10,"smiles_text":"CC(=O)Oc1ccccc1C(=O)O\nCC(=O)Nc1ccc(O)cc1"}'
