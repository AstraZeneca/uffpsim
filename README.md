# Ultra-fast Chemistry search
A Python library implementing algorithm for searching large chemistry database with ultra-fast performance.

# Installation
The simplest way to install it, is through conda environment. Following will create a new conda environment in 
current working directory and will install all the requirements from `dev-environment.yaml` file.

    # Create development conda environment
    conda env create --prefix ./venv --file dev-environment.yaml

    # activate the environment
    conda activate ./venv

    # install uffpsim
    pip install -v .

## Install with `-march=native` support
This will compile and build `uffpsim` that will be highly optimized for the current CPU. 

    UFFPSIM_NATIVE=1 pip install -v .

## Install with AVX512-VPOPCNTDQ support
The [AVX512-VPOPCNTDQ](https://en.wikipedia.org/wiki/AVX-512#VPOPCNTDQ_and_BITALG) instruction set could speed-up similarity calculation by more than 50% on newest Intel CPUs.

To check whether CPU support AVX512-VPOPCNTDQ, following command in linux should show `avx512_vpopcntdq` in the list.

    lscpu | grep Flags

If `avx512_vpopcntdq` is not in the above list, it means CPU does not have the capability to use this method.

This will compile and build `uffpsim` with AVX512-VPOPCNTDQ suppport. 

    UFFPSIM_AVX512=1 pip install -v .



# How to use?

## Creating database

### Creating database in serial with inner-clustering also in serial
Following snippet shows, how to create database where IDs are created automatically (serial number). `mol_id_max_chars` is number of maximum charecters in molecule-id, and it should be multiple of 8 minus one, e.g. 7, 15, 23, 31 etc.

```python
from uffpsim import create_database

create_database(input_file='chembl_33.sdf.gz', db_file='chembl_2048b.h5', fp_type='Morgan', fp_params={'fpSize': 2048, 'radius': 2}, gen_ids = True, mol_id_prop=None, mol_id_max_chars= 15, inner_clustering_threshold=0.15)

```

### Creating database in serial with inner-clustering in parallel

The number of threads in inner-clustering can be controlled by setting `OMP_NUM_THREADS` environment variable e.g. `export OMP_NUM_THREADS=4`.

```python
from uffpsim import create_database

create_database(input_file='chembl_33.sdf.gz', db_file='chembl_2048b.h5', fp_type='Morgan', fp_params={'fpSize': 2048, 'radius': 2}, gen_ids = True, mol_id_prop=None, mol_id_max_chars= 15, inner_clustering_threshold=0.15, cluster_parallel=True)

```

### Creating database in parallel with inner-clustering in parallel

`create_database_parallel` accepts `workers`, which can be used to set the number of threads for creating database in parallel.

The number of threads in inner-clustering can be controlled by setting `OMP_NUM_THREADS` environment variable e.g. `export OMP_NUM_THREADS=4`.

Note: creating database in parallel and inner-clustering in parallel uses two different mechanism, therefore requires two different parameters to control it.

```python
from uffpsim import create_database_parallel

create_database_parallel(input_file='chembl_33.sdf.gz', db_file='chembl_2048b.h5', workers=4, fp_type='Morgan', fp_params={'fpSize': 2048, 'radius': 2}, gen_ids = True, mol_id_prop=None, mol_id_max_chars= 15, inner_clustering_threshold=0.15, cluster_parallel=True)

```

### Redoing inner-clustering

Inner-clustering can be re-performed with a different clustering threshold value. It is good to check performance with different threshold values for new database.
This step will write new clustering data in same input file.

```python
from uffpsim import redo_inner_clustering

redo_inner_clustering("chembl_2048b.h5", 0.1, cluster_parallel=True)

```

## Searching database

### Searching database sequentially

Each time only one smiles can be used as input and it returns the result.

```python
from uffpsim import UFFPSimSearchEngine

search_engine = UFFPSimSearchEngine("chembl_2048b.h5")

# print few parameters from database
print(search_engine.fp_store.fp_params_json)
print(search_engine.fp_store.mol_id_max_chars)
print(search_engine.fp_store.fp_bits_size)
print(search_engine.fp_store.inner_clustering_threshold)

# similarity threshold of 0.8 and return only upto one hit
result_1 = search_engine.search("Cc1cc(-n2ncc(=O)[nH]c2=O)ccc1C(=O)c1ccc(C#N)cc1", 0.8, limits=1)
print(result_1)

# similarity threshold of 0.6 and return only upto 10 hits
result_2 = search_engine.search("Cc1cc(-n2ncc(=O)[nH]c2=O)ccc1C(=O)c1ccc(C#N)cc1", 0.6, limits=10)
print(result_2)

```

### Searching database in batch mode

Searching in batch mode could potentially speed-up the search by upto two times depending on the nature of queries.

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
result_1 = search_engine.batch_search(smiles_list, 0.8, limits=1)
print(result_1)

# similarity threshold of 0.6 and return only upto 10 hits
result_2 = search_engine.batch_search(smiles_list, 0.6, limits=10)
print(result_2)

```


## Development Setup

    conda env create --prefix ./venv --file dev-environment.yaml # Create development conda environment
    conda activate ./venv
    pip install -ve .
