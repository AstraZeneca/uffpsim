"""Search every query with each method and compare the hits with an exhaustive exact search.

The reference computes every Tanimoto similarity as an exact fraction (shared bits / union bits),
so it involves no floating-point rounding. For each query, threshold and k, a method's result is
  identical          its hits have exactly the reference similarity values
                     (for top-k, compounds tied at the k-th value are interchangeable);
  threshold_ties     it differs only by compounds whose similarity equals the threshold exactly;
  different          anything else (missing hits, extra compounds or duplicates).
"""
import csv
import sys
from fractions import Fraction
from multiprocessing import Pool

import numpy as np
from FPSim2 import FPSim2Engine
from uffpsim import UFFPSimSearchEngine

from common import database_dir, fingerprint_words, fp_params, load_config, read_query_smiles

FIELDS = ["query", "threshold", "k", "method", "n_expected", "n_returned", "outcome", "max_score_error"]
STATE = {}


def init_worker(config, database):
    out_dir = database_dir(config, database)
    words = np.load(out_dir / "reference_words.npy", mmap_mode="r")
    STATE.update(
        config=config,
        params=fp_params(config, database),
        ids=np.load(out_dir / "reference_ids.npy"),
        words=words,
        popcounts=np.bitwise_count(words).sum(axis=1, dtype=np.int64),
        uffpsim={mode: UFFPSimSearchEngine(str(out_dir / "uffpsim.h5"), mode) for mode in config["uffpsim"]["modes"]},
        fpsim2=FPSim2Engine(str(out_dir / "fpsim2.h5")),
    )


def exhaustive_hits(query_words, threshold):
    common = np.bitwise_count(STATE["words"] & query_words).sum(axis=1, dtype=np.int64)
    union = STATE["popcounts"] + np.bitwise_count(query_words).sum(dtype=np.int64) - common
    t = Fraction(str(threshold))
    hit = common * t.denominator >= union * t.numerator
    return {int(i): Fraction(int(c), int(u)) for i, c, u in zip(STATE["ids"][hit], common[hit], union[hit])}


def search_uffpsim(mode, smiles, threshold, k):
    limit = len(STATE["ids"]) if k == "all" else k
    return [(int(mol_id), score) for mol_id, score in STATE["uffpsim"][mode].search(smiles, threshold, limit)]


def search_fpsim2(smiles, threshold, k):
    engine = STATE["fpsim2"]
    hits = engine.similarity(smiles, threshold) if k == "all" else engine.top_k(smiles, k, threshold)
    return [(int(mol_id), float(score)) for mol_id, score in hits]


def methods():
    for mode in STATE["config"]["uffpsim"]["modes"]:
        yield f"UFFPSim {mode}", lambda s, t, k, mode=mode: search_uffpsim(mode, s, t, k)
    yield "FPSim2", search_fpsim2


def classify(returned, reference, expected, threshold):
    ids = [mol_id for mol_id, _ in returned]
    if len(set(ids)) != len(ids) or any(mol_id not in reference for mol_id in ids):
        return "different"
    got = sorted((reference[mol_id] for mol_id in ids), reverse=True)
    if got == expected:
        return "identical"
    above = lambda scores: [s for s in scores if s != threshold]
    return "threshold_ties" if above(got) == above(expected) else "different"


def compare_query(item):
    index, smiles = item
    config = STATE["config"]
    query_words = fingerprint_words(smiles, config["fingerprint"]["type"], STATE["params"])
    if query_words is None:
        raise ValueError(f"RDKit cannot parse query {index}: {smiles}")
    searches = config["searches"]
    candidates = exhaustive_hits(query_words, min(searches["thresholds"]))
    rows = []
    for threshold in searches["thresholds"]:
        t = Fraction(str(threshold))
        reference = {mol_id: s for mol_id, s in candidates.items() if s >= t}
        ranked = sorted(reference.values(), reverse=True)
        for k in ["all", *searches["top_k"]]:
            expected = ranked if k == "all" else ranked[:k]
            for name, search in methods():
                returned = search(smiles, threshold, k)
                errors = [abs(score - float(reference[i])) for i, score in returned if i in reference]
                rows.append({
                    "query": index, "threshold": threshold, "k": k, "method": name,
                    "n_expected": len(expected), "n_returned": len(returned),
                    "outcome": classify(returned, reference, expected, t),
                    "max_score_error": max(errors, default=0.0),
                })
    return rows


def compare_database(config, database, queries):
    out_path = database_dir(config, database) / "per_query.csv"
    with Pool(config["workers"], init_worker, (config, database)) as pool, open(out_path, "w", newline="") as f:
        writer = csv.DictWriter(f, FIELDS)
        writer.writeheader()
        for n, rows in enumerate(pool.imap_unordered(compare_query, enumerate(queries), chunksize=8), 1):
            writer.writerows(rows)
            if n % 1000 == 0:
                print(f"b{database['fp_size']}: {n:,}/{len(queries):,} queries", flush=True)


def main():
    config = load_config(sys.argv[1])
    queries = read_query_smiles(config)
    for database in config["databases"]:
        compare_database(config, database, queries)


if __name__ == "__main__":
    main()
