# Ultrafast Fingerprint Similarity (`UFFPSim`)

**UFFPSim** is a high-performance library for exact chemical fingerprint similarity search over large molecular databases. It extends the BitBound algorithm with a second pruning stage based on clustered fingerprints within each popcount bin, reducing the number of exact Tanimoto comparisons required. The library supports both in-memory search for high-throughput screening and disk-based search for databases that exceed available RAM, enabling exact queries on databases with up to a billion compounds.

# Installation
The simplest way to install is:

    UFFPSIM_NATIVE=1 pip install -v .

If the installation fails due to missing dependencies or shared libraries, you can alternatively use a conda environment. The following will create a new conda environment in the current working directory and install all requirements from the `dev-environment.yaml` file.

    # Create development conda environment
    conda env create --prefix ./venv --file dev-environment.yaml

    # activate the environment
    conda activate ./venv

    # install uffpsim
    UFFPSIM_NATIVE=1 pip install -v .

## Install with `-march=native` support
This will compile and build `uffpsim` that will be highly optimized for the current CPU. 

    UFFPSIM_NATIVE=1 pip install -v .

## Install with AVX512-VPOPCNTDQ support
The [AVX512-VPOPCNTDQ](https://en.wikipedia.org/wiki/AVX-512#VPOPCNTDQ_and_BITALG) instruction set could speed-up similarity calculation by more than 50% on newest Intel CPUs.

To check whether CPU support AVX512-VPOPCNTDQ, following command in linux should show `avx512_vpopcntdq` in the list.

    lscpu | grep Flags | grep avx

If `avx512_vpopcntdq` is not in the above list, it means CPU does not have the capability to use this method.

This will compile and build `uffpsim` with AVX512-VPOPCNTDQ support. 

    UFFPSIM_AVX512=1 pip install -v .

# How to use?

The package can be used as a Python library or as the command.

## Command line interface (CLI)



| Command                    | Description                                                                 |
|---------------------------|-----------------------------------------------------------------------------|
| ``create-database``       | Create a new uffpsim HDF5 database from a molecular input file.             |
| ``redo-clustering``       | Redo inner clustering for an existing database with a new similarity threshold. |
| ``build-mol-id-index-table`` | Build serialized MolIdIndexTable for an existing database.              |
| ``search``                | Search one or more SMILES against the database and write hits to a CSV file.   |
| ``launch-web-app``        | Launch the uffpsim interactive web application for visual similarity searching. |

**For more details about the CLIs, visit the documentation site.**

## Creating database

```bash
# Download ChEMBL 33 to use as an example
curl https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases/chembl_33/chembl_33.sdf.gz -o ./chembl_33.sdf.gz
ls -artl chembl_33.sdf.gz
# -rw-rw-r-- 1 user user 770312723 Jun 11 12:28 chembl_33.sdf.gz
```

### Serial ID creation
IDs are created automatically as string serial numbers starting from `'1'`. `mol_id_max_chars` sets the maximum number of characters in a molecule ID and must be a multiple of 8 minus one (e.g. 7, 15, 23, 31). This value is fixed at database creation time and cannot be changed afterwards.

Molecule IDs are stored inline in the same fixed-size array as the fingerprint bits, so every entry in the database reserves exactly `mol_id_max_chars + 1` bytes for the ID regardless of actual ID length. Set it to comfortably fit your longest molecule ID — for example, ChEMBL IDs (`CHEMBL1234567`, 13 chars) fit within `mol_id_max_chars=15`.

> **Warning:** IDs longer than `mol_id_max_chars` are silently truncated, which can cause incorrect lookups. Choose a value large enough for your dataset upfront.

**Memory trade-off:** each additional 8 characters of headroom costs 8 bytes per molecule. As a ballpark example, for a 2M-molecule database, increasing from `mol_id_max_chars=15` to `mol_id_max_chars=31` adds ~16 MB to the database size.

### Creating database (serial fingerprints, serial clustering)

Fingerprints are calculated one molecule at a time and clustering runs in a single thread. This is the simplest option and uses the least memory.

```python
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
```

### Creating database (serial fingerprints, parallel clustering)

Fingerprint calculation is still serial, but the inner-clustering step uses multiple OpenMP threads within the same process. Set `OMP_NUM_THREADS` to control the thread count, e.g. `export OMP_NUM_THREADS=4`.

```python
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
```

### Creating database (parallel fingerprints, parallel clustering)

Fingerprint calculation is distributed across `workers` independent Python processes using `ProcessPoolExecutor`, which bypasses the GIL and speeds up the RDKit computation step. Writing to the HDF5 file and inner-clustering remain serial and parallel (via OpenMP) respectively. Note that spawning multiple processes requires more memory than the serial approach.

Set `OMP_NUM_THREADS` to control clustering threads independently of `workers`.

```python
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
```

### Redoing inner-clustering

Inner-clustering can be re-performed with a different threshold. It is useful to benchmark different threshold values for a new database. The clustering data is written back to the same file.

```python
from uffpsim import redo_inner_clustering

redo_inner_clustering("chembl_2048b.h5", 0.1, cluster_parallel=True)
```

## Searching database

### Searching database sequentially

Each time only one SMILES string can be used as input.

```python
from uffpsim import UFFPSimSearchEngine

search_engine = UFFPSimSearchEngine("chembl_2048b.h5")

# print few parameters from database
print(search_engine.fp_store.fp_params_json)
print(search_engine.fp_store.mol_id_max_chars)
print(search_engine.fp_store.fp_bits_size)
print(search_engine.fp_store.inner_clustering_threshold)

# similarity threshold of 0.8 and return only up to one hit
result_1 = search_engine.search("Cc1cc(-n2ncc(=O)[nH]c2=O)ccc1C(=O)c1ccc(C#N)cc1", 0.8, limit_by=1)
print(result_1)

# similarity threshold of 0.6 and return only up to 10 hits
result_2 = search_engine.search("Cc1cc(-n2ncc(=O)[nH]c2=O)ccc1C(=O)c1ccc(C#N)cc1", 0.6, limit_by=10)
print(result_2)
```

```python
# Example of result from search
>>> search_engine.search("Cc1cc(-n2ncc(=O)[nH]c2=O)ccc1C(=O)c1ccc(C#N)cc1", 0.8, limit_by=1)
[('2', 1.0)]  # '2' matches the second SMILES string in `uffpsim/tests/data/10mols.smi`
```

### Searching database in batch mode

Searching in batch mode could potentially speed-up the search by up to two times depending on the nature of queries.

```python
from uffpsim import UFFPSimSearchEngine

search_engine = UFFPSimSearchEngine("chembl_2048b.h5")

# list of smiles to be searched
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

# similarity threshold of 0.8 and return only upto one hit
result_1 = search_engine.batch_search(smiles_list, 0.8, limit_by=1)
print(result_1)

# similarity threshold of 0.6 and return only upto 10 hits
result_2 = search_engine.batch_search(smiles_list, 0.6, limit_by=10)
print(result_2)
```

## Development Setup
    conda env create --prefix ./venv --file dev-environment.yaml # Create development conda environment
    conda activate ./venv
    pip install -ve .
