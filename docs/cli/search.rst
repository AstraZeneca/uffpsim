search
======

Search one or more SMILES strings against a database and write hits to a CSV output file.

.. code-block:: bash

    uffpsim search -d <db_file> -o <output_csv> [options]

Options
-------

.. program:: uffpsim search

.. option:: -d, --db-file <path>

    (Required) HDF5 database path to search against.

.. option:: -q, --smiles <SMILES>

    Single SMILES query string from the command line.

.. option:: -i, --smiles-file <path>

    Path to an input file containing one SMILES per line. Lines starting with ``#`` are treated as comments and skipped.

.. option:: -o, --output-csv <path>

    (Required) Output CSV file path for search results.

.. option:: -m, --mode {memory,disk}

    Search engine load mode. ``memory`` loads the database into RAM; ``disk`` accesses data from the HDF5 file directly. Default: memory.

.. option:: -t, --threshold <float>

    Similarity threshold (Tanimoto coefficient). Default: 0.6.

.. option:: -k, --limit <int>

    Maximum number of hits to return per query. Default: 10.

.. option:: -s, --include-hit-smiles

    Include hit SMILES in the output CSV alongside compound IDs and scores.

Output CSV Format
-----------------

The output CSV columns are:

+-----------------+-----------------------------------------------------------+
| Column          | Description                                               |
+=================+===========================================================+
| query_index     | Sequential index of the query (1-based)                   |
+-----------------+-----------------------------------------------------------+
| query_smiles    | The input SMILES string                                   |
+-----------------+-----------------------------------------------------------+
| hit_rank        | Rank among hits for this query (1 = highest score)        |
+-----------------+-----------------------------------------------------------+
| compound_id     | Database compound identifier                              |
+-----------------+-----------------------------------------------------------+
| score           | Tanimoto similarity score                                 |
+-----------------+-----------------------------------------------------------+
| hit_smiles      | SMILES of the hit (only if ``--include-hit-smiles`` set)  |
+-----------------+-----------------------------------------------------------+
| error           | Error message if fingerprint generation failed            |
+-----------------+-----------------------------------------------------------+

Examples
--------

Single query with default threshold:

.. code-block:: bash

    uffpsim search -d chembl.h5 -q "Cc1cc(-n2ncc(=O)[nH]c2=O)ccc1C(=O)c1ccc(C#N)cc1" \
        -o results.csv --limit 5

Batch query from a file, with hit SMILES and custom threshold:

.. code-block:: bash

    uffpsim search -d chembl.h5 -i queries.smi -o batch_results.csv \
        --threshold 0.5 --limit 20 --include-hit-smiles
