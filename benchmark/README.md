# Benchmarking UFFPSim Performance

This directory contains the scripts and data used to reproduce the benchmarking experiments presented in the UFFPSim manuscript.

The benchmarks cover:

- **Query selection** – construction of the ZINC query set used throughout the experiments.
- **ChEMBL-35** – fingerprint-size and cluster-threshold optimisation, and comparison with FPSim2.
- **PubChem** – comparison of the 'static-threshold' and 'stepped' search strategies, in-memory and disk-based implementations, and FPSim2.
- **ZINC20** – scalability evaluation of the disk-based implementation on a nearly one-billion-compound database.

## Installation

The benchmarks use the `sim-ops` branch of UFFPSim:

```bash
conda env create --prefix ./venv --file dev-environment.yaml
conda activate ./venv
git checkout sim-ops
UFFPSIM_NATIVE=1 pip install .
```

The `sim-ops` branch additionally reports the number of fingerprints scanned during a search, which is used for the search-strategy analysis.

## Query selection

The benchmark query set is provided in:

```text
zinc_query_smiles.txt
```

The query-selection procedure is described in the manuscript.

## ChEMBL-35

Scripts in this directory reproduce the cluster membership threshold optimisation and the UFFPSim/FPSim2 comparison across fingerprint sizes and top-*k* search settings.

UFFPSim/FPSim2 comparisons use the **'static-threshold' in-memory implementation**.

## PubChem

Scripts reproduce the PubChem benchmark comparing:

- **'static-threshold' in-memory**
- **'stepped' in-memory**
- **disk-based 'stepped'**

UFFPSim configurations are also compared with FPSim2. The search-strategy analysis uses the number of fingerprints scanned by the two in-memory strategies.

## ZINC20

Scripts reproduce the scalability benchmark using the disk-based UFFPSim implementation on the nearly one-billion-compound ZINC20 subset.

FPSim2 was not benchmarked at this scale because the complete fingerprint database could not reside in memory on the benchmark hardware.

## Compilation options

UFFPSim supports architecture-specific compilation options, including:

```bash
UFFPSIM_NATIVE=1
UFFPSIM_AVX512=1
```

The AVX512-VPOPCNTDQ configuration was not used for the benchmark results reported in the manuscript.

Further details of the benchmark procedures, datasets, parameters, and hardware are provided in the UFFPSim manuscript and Supporting Information.
