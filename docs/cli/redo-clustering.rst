redo-clustering
===============

Redo inner clustering for an existing database with a new similarity threshold.

.. code-block:: bash

    uffpsim redo-clustering -d <db_file> -t <threshold> [options]

Options
-------

.. option:: -d, --db-file <path>

    (Required) HDF5 database path to re-cluster.

.. option:: -t, --threshold <float>

    (Required) New inner clustering threshold. Molecules above this Tanimoto similarity will be grouped together.

.. option:: -c, --cluster-mode {memory,disk}

    Clustering mode. ``memory`` stores clustering data in RAM only; ``disk`` writes it to the database file. Default: memory.

.. option:: -P, --cluster-parallel

    Enable OpenMP-based clustering parallelism. Controls thread count via ``OMP_NUM_THREADS``.

Example
-------

Redo clustering with a lower threshold for finer granularity:

.. code-block:: bash

    uffpsim redo-clustering -d chembl.h5 -t 0.1 --cluster-parallel
