Demo Web-App
============

This page provides guides for running uffpsim's built-in web application locally in two demo scenarios.

Prerequisites
-------------

Before starting, install uffpsim and its dependencies:

.. code-block:: bash

    UFFPSIM_NATIVE=1 pip install -v ".[web]"

Demo 1: Quick Search with a Public Dataset Snapshot
-----------------------------------------------------

Step 1: Download a Database File
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For demonstration purposes, you can create a small database from the test data included in uffpsim:

.. code-block:: python

    # From within the uffpsim repository:
    from uffpsim import create_database

    create_database(
        input_file='tests/data/10mols.sdf',
        db_file='demo_db.h5',
        fp_type='Morgan',
        fp_params={'fpSize': 2048, 'radius': 2},
        gen_ids=True,
        inner_clustering_threshold=0.15,
    )

Step 2: Launch the Web App
~~~~~~~~~~~~~~~~~~~~~~~~~~

Run the following command from the uffpsim repository root:

.. code-block:: bash

    python -m uffpsim launch-web-app -d demo_db.h5

Or via the installed entry point:

.. code-block:: bash

    uffpsim launch-web-app -d demo_db.h5

Step 3: Access the Web Application
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Open a browser and navigate to:

    http://127.0.0.1:5000/

You can now enter SMILES strings in the search form and view similarity results with molecular structure images.

Demo 2: Full Search with a Large ChEMBL Database
--------------------------------------------------

Step 1: Create or Download a Large Database
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For this demo, create a database from a larger molecular dataset such as ChEMBL:

.. code-block:: bash

    # Download ChEMBL 33 (example)
    curl https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases/chembl_33/chembl_33.sdf.gz -o ./chembl_33.sdf.gz

.. code-block:: python

    from uffpsim import create_database_parallel

    create_database_parallel(
        input_file='chembl_33.sdf.gz',
        db_file='chembl_2048b.h5',
        workers=4,
        fp_type='Morgan',
        fp_params={'fpSize': 2048, 'radius': 2},
        gen_ids=True,
        mol_id_max_chars=15,
        inner_clustering_threshold=0.15,
        cluster_parallel=True,
    )

Step 2: Launch the Web App with Full Results Mode
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Use ``--results full`` to return compound IDs, SMILES strings, and molecular images for each hit.
Build the molecule ID index first so SMILES lookup is available:

.. code-block:: bash

    uffpsim build-mol-id-index-table -d chembl_2048b.h5

    python -m uffpsim launch-web-app \
        -d chembl_2048b.h5 \
        --results full \
        --port 5000

Or with memory mode explicitly set:

.. code-block:: bash

    uffpsim launch-web-app -d chembl_2048b.h5 -m memory -r full -p 5000

Step 3: Access the Web Application
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Navigate to:

    http://127.0.0.1:5000/

The web app will display full hit details including molecular structures for each matching compound.

Web Application Options
-----------------------

.. list-table::
   :widths: 25 15 60
   :header-rows: 1

   * - Option
     - Default
     - Description
   * - ``-d, --db-file``
     - (required)
     - Path to the HDF5 database file to preload at startup.
   * - ``-m, --mode``
     - ``memory``
     - Search engine load mode: ``memory`` loads data into RAM; ``disk`` keeps data in the database file.
   * - ``-r, --results``
     - ``ids_only``
     - Result detail level: ``ids_only`` returns compound IDs and scores only; ``full`` also returns SMILES strings and RDKit images.
   * - ``-H, --host``
     - ``127.0.0.1``
     - Host interface to bind the Flask server to.
   * - ``-p, --port``
     - ``5000``
     - TCP port to bind the Flask server to.
   * - ``-D, --debug``
     - Disabled
     - Enable Flask debug mode for development.
