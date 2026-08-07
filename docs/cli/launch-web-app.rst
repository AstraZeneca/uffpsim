launch-web-app
==============

Launch the uffpsim interactive web application for visual similarity searching.

.. code-block:: bash

    uffpsim launch-web-app -d <db_file> [options]

Options
-------

.. option:: -d, --db-file <path>

    (Required) HDF5 database path to preload at startup.

.. option:: -m, --mode {memory,disk}

    Search engine load mode. ``memory`` loads data into RAM for faster queries; ``disk`` accesses the database file directly. Default: memory.

.. option:: -r, --results {ids_only,full}

    Result detail level. ``ids_only`` returns compound IDs and scores; ``full`` also fetches SMILES strings and molecular images (requires MolIdIndexTable to be built first). Default: ids_only.

.. option:: -H, --host <ip>

    Host interface to bind. Default: 127.0.0.1.

.. option:: -p, --port <int>

    TCP port to bind. Default: 5000.

.. option:: -D, --debug

    Enable Flask debug mode for development.

Example
-------

Launch with default settings:

.. code-block:: bash

    uffpsim launch-web-app -d chembl.h5

Launch with full results and custom port:

.. code-block:: bash

    uffpsim launch-web-app -d chembl.h5 --results full -p 8080

Access the web app at http://127.0.0.1:5000/ (or ``http://<host>:<port>/``).
