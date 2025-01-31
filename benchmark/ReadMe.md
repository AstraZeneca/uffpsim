## Benchmarking UFFPSim performance

### Slow vs Fast Search Query
The algorithm implemented in [FPSim2](https://github.com/chembl/FPSim2) is fastest if hits are found with very high 
similarity, e.g. `1.0`. However, the search time increases with decrease in maximum similarity of hits due to the
ever widening of search space. Therefore, the performance of search query also depends on how much a compound is most
similar in the target database.

The counting of 1's in fingerprint is extremely optimized and its performance depends on number of 1's in the fingerprint.
Therefore, query with large number of 1's in fingerprint might be slower than that of fingerprint with fewer 1's.

### Preparation of queries and target dataset
As discussed above in Slow vs Fast Search Query, the benchmark queries should have following two properties

* It should have both slow and fast queries - low to 1.0 similarity score with target database
* 1's count in fingerprints should span a wide range

Based on the above two considerations, we prepared the benchmark queries as follows.
A list of queries was prepared with diverse similarity to target database. Here, we selected queries
from ZINC database and Chembl-33 as target database. The ZINC database was screened against Chembl-33
with 2048 bits of Morgan fingerprint and radius of 2 $\AA$. The final queries from ZINC database are selected such
that their maximum similarity in Chembl-33 ranges from 0.5 to 1.0. Additionally, the one's count in
queries fingerprint also ranges from ~5 to ~100. Therefore, our aim was to diversify similarity scores and
one's count among queries.

The selected queries are listed in `zinc_query_smiles.txt` with ZINC-ID and highest similarity score in Chembl-33.

### Benchmarking

#### With default compilation
When default compilation is used to build UFFPSim, it uses SSE4.2 instruction set is used during compilation.
Additionally `_builtin_popcount` is used to count number of one's bits in the fingerprints.

#### With `UFFPSIM_NATIVE=1` environment variable
When this variable is used, it uses `-march=native` flag and build the most optimized version for the that machine.
This version also uses `_builtin_popcount` to count number of one's bits in the fingerprints.
There might be possibility that when compiled and built code is used on another machine, its performance could
degrade.

#### With `UFFPSIM_AVX512=1` environment variable
When this variable is used, it uses [AVX512-VPOPCNTDQ](https://en.wikipedia.org/wiki/AVX-512#VPOPCNTDQ_and_BITALG) instruction
set to count number of one's bits in the fingerprints together with `AND` operation.

To check whether CPU support AVX512-VPOPCNTDQ, following command in linux should show `avx512_vpopcntdq` in the list.

    lscpu | grep Flags

If `avx512_vpopcntdq` is not in the above list, it means CPU does not have the capability to use this method.

### Batch vs Sequential

The speed-up factor for batch-mode was calculated against the sequential mode. For this, randomly batch of queries 
were selected and search-time was computed in batch-search and sequential-search mode. Subsequently, speed-up
factor for each batch of queries were obtained.

The speed-up depends on the nature of queries in the given batch. Particularly, it depends of the count of 1's in the 
fingerprints. If queries have very similar 1's count, performance gain is expected to be larger. However, as 1's
count become more diverse, performance might be worse than the sequential-search mode.
