# 01 – Exactness: UFFPSim vs FPSim2

Do UFFPSim's inner clustering and cluster screening change which compounds a search returns?

UFFPSim (all three search modes) and FPSim2 are run on the same queries and the same fingerprints,
and every result is compared with an exhaustive search that computes each Tanimoto similarity as an
exact fraction (shared bits / union bits), without floating-point rounding. A method that agrees
with this reference for every query returns exactly the compounds an exhaustive search would; if
UFFPSim and FPSim2 both agree with it, they agree with each other.

## Setup

| | |
|---|---|
| Database | ChEMBL 35 ([metadata/chembl_35.yaml](metadata/chembl_35.yaml)) |
| Queries | the 22,581 ZINC benchmark queries of the manuscript ([metadata/zinc_queries.yaml](metadata/zinc_queries.yaml)) |
| Fingerprint | Morgan, radius 2; sizes and inner-clustering thresholds as in the manuscript ([config.yaml](config.yaml)) |
| UFFPSim modes | `memory` (static threshold), `memory-stepped`, `disk` |
| Searches | all hits at thresholds 0.5, 0.6, 0.7; top-1 and top-10 at each threshold |

UFFPSim, FPSim2 and the reference are built from the same SMILES file, and the build fails unless the
FPSim2 fingerprints are bit-for-bit identical to the reference fingerprints, so the comparison tests
the search algorithms only.

## Outcomes

For each query, threshold and *k*, a result is

- **identical** – its hits have exactly the reference similarity values. For top-*k*, compounds tied
  at the *k*-th value are interchangeable, so only the similarity values are compared;
- **threshold_ties** – it differs from the reference only by compounds whose similarity equals the
  threshold exactly;
- **different** – anything else (a missing hit, an extra compound or a duplicate).

## Running

From the repository root:

```bash
conda env create --prefix ./venv-exp01 --file experiments/01-exactness_comparison_fpsim2/environment.yaml
conda activate ./venv-exp01
UFFPSIM_NATIVE=1 pip install .

cd experiments/01-exactness_comparison_fpsim2
curl -L -o data/chembl_35_chemreps.txt.gz "$(python -c 'import yaml; print(yaml.safe_load(open("metadata/chembl_35.yaml"))["url"])')"
python build_databases.py config.yaml   # SMILES file, reference fingerprints, UFFPSim and FPSim2 databases
python compare_searches.py config.yaml  # data/b<size>/per_query.csv
python summarise.py config.yaml         # results/summary.md
```

`data/` is not tracked. The summary table is in [results/summary.md](results/summary.md).
