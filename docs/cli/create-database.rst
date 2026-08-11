create-database
===============

Create a new uffpsim HDF5 database from a molecular input file or from all matching files in a directory.

.. code-block:: bash

    uffpsim create-database -i <input_file> -d <db_file> -f <fp_type> [options]
    uffpsim create-database -D <input_dir> -d <db_file> -f <fp_type> [options]

Options
-------

.. program:: uffpsim create-database

.. option:: -i, --input-file <path>

    Input molecule file path. Supported formats: SDF, SDF.GZ, SMILES (SMI), InChI.

.. option:: -D, --input-dir <path>

    Input directory containing molecular files. When provided, all matching files in the directory are processed. Requires ``-s/--suffix`` to choose the file type.

.. option:: -d, --db-file <path>

    (Required) Output HDF5 database path.

.. option:: -f, --fp-type <type>

    (Required) Fingerprint type. One of: Morgan, RDKit, AtomPair, TopologicalTorsion, MACCSKeys, Avalon, RDKPatternFingerprint.

.. option:: -p, --fp-params <JSON>

    Fingerprint parameters as a JSON object. Example: ``{"fpSize": 2048, "radius": 2}``. Defaults to fingerprint type defaults.

.. option:: -w, --workers <int>

    Number of workers. Use >1 for parallel processing mode (spawns multiple Python processes). Default: 1.

.. option:: -g, --gen-ids / --no-gen-ids

    Generate IDs if missing in input. If the input file already contains IDs, set ``--no-gen-ids``. Default: ``--gen-ids``.

.. option:: -m, --mol-id-prop <name>

    Property name to use as molecule ID from the input file (e.g., a field in an SDF file). Defaults to None (auto-generated serial IDs).

.. option:: -l, --mol-id-max-chars <int>

    Maximum molecule ID length. Must be a multiple of 8 minus one (7, 15, 23, 31...). Default: 15.

.. option:: -I, --info <text or JSON>

    Additional database info stored in the HDF5 file as a JSON string or plain text.

.. option:: -t, --inner-clustering-threshold <float>

    Threshold for inner clustering. Molecules above this Tanimoto similarity are grouped in pop-count bins. Default: 0.2.

.. option:: -c, --cluster-mode {memory,disk,cuda}

    Clustering mode. ``memory`` stores clustering data in RAM only; ``disk`` writes it to the database file; ``cuda`` uses GPU acceleration. Default: memory.

.. option:: -P, --cluster-parallel

    Enable OpenMP-based clustering parallelism. Controls thread count via ``OMP_NUM_THREADS``.

.. option:: -s, --suffix <suffix>

    File suffix to scan when using ``--input-dir``. Supported values: ``sdf``, ``sdf.gz``, ``smi``. Default: ``sdf.gz``.

Example
-------

Create a database from an SDF file with Morgan fingerprints:

.. code-block:: bash

    uffpsim create-database \
        -i chembl_33.sdf.gz \
        -d chembl.h5 \
        -f Morgan \
        --fp-params '{"fpSize": 2048, "radius": 2}' \
        --gen-ids \
        --mol-id-max-chars 15 \
        --inner-clustering-threshold 0.15

Create a database using parallel fingerprint processing with 4 workers:

.. code-block:: bash

    uffpsim create-database \
        -i chembl_33.sdf.gz \
        -d chembl.h5 \
        -f Morgan \
        --fp-params '{"fpSize": 2048, "radius": 2}' \
        -w 4 \
        --cluster-parallel

Create a database from all SMILES files in a directory:

.. code-block:: bash

    uffpsim create-database \
        -D /path/to/molecules \
        -d chembl.h5 \
        -f Morgan \
        --suffix smi \
        --gen-ids
