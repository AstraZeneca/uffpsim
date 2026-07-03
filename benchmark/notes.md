# Benchmarking UFFPSim Performance

## Slow vs Fast Search Query

The algorithm implemented in [FPSim2](https://github.com/chembl/FPSim2) is fastest when hits are found with very high similarity, for example `1.0`. As the maximum similarity of the hits decreases, the search space widens, resulting in longer search times. Therefore, search performance depends on the similarity of a query molecule to its closest match in the target database.

The counting of set bits in a fingerprint is highly optimised, but its performance still depends on the number of set bits. Queries with a larger number of set bits may therefore take longer to process than those with fewer set bits.

## Preparation of Queries and Target Dataset

The benchmark query set was designed to satisfy two requirements:

- include both fast and slow queries by spanning a range of maximum similarities to the target database, from low similarity to `1.0`
- include queries with a wide range of fingerprint populations

To achieve this, query molecules were selected from the ZINC database, with ChEMBL-33 used as the target database. The ZINC database was screened against ChEMBL-33 using 2048-bit Morgan fingerprints with a radius of 2. The final query set was selected such that the maximum similarity of each query to ChEMBL-33 ranged from approximately `0.5` to `1.0`. Additionally, the number of set bits in the query fingerprints ranged from approximately `5` to `100`.

The selected queries are provided in `zinc_query_smiles.txt`, together with their ZINC identifiers and highest similarity score in ChEMBL-33.

### Benchmarking

#### Default Compilation

By default, UFFPSim is compiled using the SSE4.2 instruction set. The implementation uses `__builtin_popcount` to count the number of set bits in fingerprints.

#### `UFFPSIM_NATIVE=1`

When the `UFFPSIM_NATIVE=1` environment variable is set, UFFPSim is compiled with the `-march=native` compiler flag, producing a binary optimized for the local machine. This build also uses `__builtin_popcount` to count the number of set bits in fingerprints.

As the generated binary is optimized for the architecture of the build machine, performance may degrade when it is executed on a different machine.

#### `UFFPSIM_AVX512=1`

When the `UFFPSIM_AVX512=1` environment variable is set, UFFPSim uses the [AVX512-VPOPCNTDQ](https://en.wikipedia.org/wiki/AVX-512#VPOPCNTDQ_and_BITALG) instruction set to perform the bitwise `AND` operation and count the number of set bits in a single step.

To check whether a CPU supports AVX512-VPOPCNTDQ, run:

```bash
lscpu | grep Flags
```

If `avx512_vpopcntdq` appears in the output, the processor supports this instruction set. Otherwise, this implementation cannot be used.

### Batch vs Sequential Search

The speed-up of batch search was measured relative to sequential search. Random batches of query molecules were selected, and the search time was measured using both batch and sequential search modes. The speed-up factor was then calculated for each batch.

The achievable speed-up depends on the composition of the query batch, particularly the number of set bits in the fingerprints. Batches containing queries with similar fingerprint populations are expected to achieve greater speed-up. As the variation in the number of set bits increases, the benefit of batch search decreases and may even become slower than sequential search.
