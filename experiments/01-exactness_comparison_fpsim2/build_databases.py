"""Build the UFFPSim, FPSim2 and exhaustive-reference databases defined in config.yaml.

All three are built from the same SMILES file, and the FPSim2 fingerprints are checked to be
identical to the reference fingerprints, so the experiment compares search algorithms only.
"""
import csv
import gzip
import sys
from functools import partial
from multiprocessing import Pool

import numpy as np
from FPSim2 import FPSim2Engine
from FPSim2.io import create_db_file
from uffpsim import create_database_parallel

from common import data_path, database_dir, fingerprint_words, fp_params, load_config


def write_database_smiles(config, path):
    dataset = config["dataset"]
    with gzip.open(data_path(config, dataset["file"]), "rt") as src, open(path, "w") as dst:
        for row in csv.DictReader(src, delimiter="\t"):
            mol_id = int(row[dataset["id_column"]].removeprefix(dataset["id_prefix"]))
            dst.write(f"{row[dataset['smiles_column']]} {mol_id}\n")


def _reference_row(line, fp_type, params):
    smiles, mol_id = line.split()
    words = fingerprint_words(smiles, fp_type, params)
    return None if words is None else (int(mol_id), words)


def write_reference(config, database, smiles_path, out_dir):
    row = partial(_reference_row, fp_type=config["fingerprint"]["type"], params=fp_params(config, database))
    with open(smiles_path) as f, Pool(config["workers"]) as pool:
        rows = [r for r in pool.imap(row, f, chunksize=5000) if r is not None]
    np.save(out_dir / "reference_ids.npy", np.array([mol_id for mol_id, _ in rows], dtype=np.int64))
    np.save(out_dir / "reference_words.npy", np.stack([words for _, words in rows]))


def build_uffpsim(config, database, smiles_path, out_dir):
    create_database_parallel(
        str(smiles_path), str(out_dir / "uffpsim.h5"), config["fingerprint"]["type"],
        workers=config["workers"], fp_params=fp_params(config, database), gen_ids=False,
        mol_id_prop=None, mol_id_max_chars=config["uffpsim"]["mol_id_max_chars"],
        inner_clustering_threshold=database["inner_clustering_threshold"],
        cluster_mode="memory", cluster_parallel=True,
    )


def build_fpsim2(config, database, smiles_path, out_dir):
    create_db_file(
        str(smiles_path), str(out_dir / "fpsim2.h5"), mol_format="smiles",
        fp_type=config["fingerprint"]["type"], fp_params=fp_params(config, database),
    )


def check_fpsim2_matches_reference(out_dir):
    fps = FPSim2Engine(str(out_dir / "fpsim2.h5")).fps
    fps = fps[np.argsort(fps[:, 0], kind="stable")]
    ids = np.load(out_dir / "reference_ids.npy")
    order = np.argsort(ids, kind="stable")
    words = np.load(out_dir / "reference_words.npy")[order]
    if not (np.array_equal(fps[:, 0].astype(np.int64), ids[order]) and np.array_equal(fps[:, 1:-1], words)):
        raise ValueError(f"FPSim2 fingerprints differ from the reference fingerprints in {out_dir}")
    print(f"{out_dir.name}: {len(ids):,} FPSim2 fingerprints identical to the reference")


def main():
    config = load_config(sys.argv[1])
    smiles_path = data_path(config, "database.smi")
    write_database_smiles(config, smiles_path)
    for database in config["databases"]:
        out_dir = database_dir(config, database)
        out_dir.mkdir(parents=True, exist_ok=True)
        write_reference(config, database, smiles_path, out_dir)
        build_uffpsim(config, database, smiles_path, out_dir)
        build_fpsim2(config, database, smiles_path, out_dir)
        check_fpsim2_matches_reference(out_dir)


if __name__ == "__main__":
    main()
