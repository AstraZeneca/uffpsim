build-mol-id-index-table
========================

Build a serialized MolIdIndexTable for an existing database to enable fast SMILES retrieval by compound ID.

.. code-block:: bash

    uffpsim build-mol-id-index-table -d <db_file>

Options
-------

.. option:: -d, --db-file <path>

    (Required) HDF5 database path for which to build the index table.

Example
-------

Build a molecule ID index for an existing ChEMBL database:

.. code-block:: bash

    uffpsim build-mol-id-index-table -d chembl.h5

After building the index, the web application's ``full`` results mode can efficiently look up SMILES strings and molecular images for each hit.
