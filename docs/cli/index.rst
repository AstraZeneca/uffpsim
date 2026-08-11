Command-Line Interface
======================

.. toctree::
   :maxdepth: 1
   :caption: Commands:
   :hidden:

   create-database
   redo-clustering
   build-mol-id-index-table
   search
   launch-web-app

Summary
-------

.. list-table::
   :widths: 30 70
   :header-rows: 1

   * - Command
     - Description
   * - :doc:`create-database <create-database>`
     - Create a new uffpsim HDF5 database from a molecular input file.
   * - :doc:`redo-clustering <redo-clustering>`
     - Redo inner clustering for an existing database with a new similarity threshold.
   * - :doc:`build-mol-id-index-table <build-mol-id-index-table>`
     - Build serialized MolIdIndexTable for an existing database.
   * - :doc:`search <search>`
     - Search one or more SMILES against the database and write hits to a CSV file.
   * - :doc:`launch-web-app <launch-web-app>`
     - Launch the uffpsim interactive web application for visual similarity searching.

Usage
-----

The CLI entry point is ``uffpsim``:

.. code-block:: bash

    uffpsim <subcommand> [options]

Type ``uffpsim --help`` to see available subcommands and their options.
