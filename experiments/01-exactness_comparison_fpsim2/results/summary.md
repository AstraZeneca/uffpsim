# Agreement with an exhaustive exact search

Each cell: queries whose hits are identical to the exhaustive reference / queries searched.
Queries with hits: queries for which the reference returns at least one compound.
Deviations, if any, are listed in brackets (see compare_searches.py for definitions).

### 2048-bit fingerprints (inner-clustering threshold 0.16)

| Search | Threshold | Queries with hits | UFFPSim memory | UFFPSim memory-stepped | UFFPSim disk | FPSim2 |
|---|---|---|---|---|---|---|
| all hits | 0.5 | 22,581 | 22,581 / 22,581 (100.00%) | 22,581 / 22,581 (100.00%) | 22,581 / 22,581 (100.00%) | 22,581 / 22,581 (100.00%) |
| top-1 | 0.5 | 22,581 | 22,581 / 22,581 (100.00%) | 22,581 / 22,581 (100.00%) | 22,581 / 22,581 (100.00%) | 22,581 / 22,581 (100.00%) |
| top-10 | 0.5 | 22,581 | 22,581 / 22,581 (100.00%) | 22,581 / 22,581 (100.00%) | 22,581 / 22,581 (100.00%) | 22,581 / 22,581 (100.00%) |
| all hits | 0.6 | 19,827 | 22,581 / 22,581 (100.00%) | 22,581 / 22,581 (100.00%) | 22,581 / 22,581 (100.00%) | 22,581 / 22,581 (100.00%) |
| top-1 | 0.6 | 19,827 | 22,581 / 22,581 (100.00%) | 22,581 / 22,581 (100.00%) | 22,581 / 22,581 (100.00%) | 22,581 / 22,581 (100.00%) |
| top-10 | 0.6 | 19,827 | 22,581 / 22,581 (100.00%) | 22,581 / 22,581 (100.00%) | 22,581 / 22,581 (100.00%) | 22,581 / 22,581 (100.00%) |
| all hits | 0.7 | 15,588 | 22,581 / 22,581 (100.00%) | 22,581 / 22,581 (100.00%) | 22,581 / 22,581 (100.00%) | 22,581 / 22,581 (100.00%) |
| top-1 | 0.7 | 15,588 | 22,581 / 22,581 (100.00%) | 22,581 / 22,581 (100.00%) | 22,581 / 22,581 (100.00%) | 22,581 / 22,581 (100.00%) |
| top-10 | 0.7 | 15,588 | 22,581 / 22,581 (100.00%) | 22,581 / 22,581 (100.00%) | 22,581 / 22,581 (100.00%) | 22,581 / 22,581 (100.00%) |

Largest difference between a reported and the exact similarity: 2.96e-08
