"""Configuration, paths and fingerprints shared by the experiment scripts."""
from functools import lru_cache
from pathlib import Path

import numpy as np
import yaml
from rdkit import Chem, RDLogger
from rdkit.Chem import rdFingerprintGenerator

EXPERIMENT_DIR = Path(__file__).resolve().parent
RDLogger.DisableLog("rdApp.*")


def load_yaml(path):
    with open(EXPERIMENT_DIR / path) as f:
        return yaml.safe_load(f)


def load_config(path):
    config = load_yaml(path)
    config["dataset"] = load_yaml(config["dataset_metadata"])
    config["query_set"] = load_yaml(config["query_metadata"])
    return config


def data_path(config, name):
    return EXPERIMENT_DIR / config["data_dir"] / name


def database_dir(config, database):
    return data_path(config, f"b{database['fp_size']}")


def fp_params(config, database):
    return {**config["fingerprint"]["params"], "fpSize": database["fp_size"]}


def read_query_smiles(config):
    query_set = config["query_set"]
    with open(EXPERIMENT_DIR / query_set["file"]) as f:
        return [line.split()[query_set["smiles_column"]] for line in f]


@lru_cache
def _generator(fp_type, params):
    return getattr(rdFingerprintGenerator, f"Get{fp_type}Generator")(**dict(params))


def fingerprint_words(smiles, fp_type, params):
    """Fingerprint as uint64 words in FPSim2's layout, or None if RDKit cannot parse the SMILES."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    bits = _generator(fp_type, tuple(sorted(params.items()))).GetFingerprintAsNumPy(mol)
    return np.packbits(bits).view(">u8").astype(np.uint64)
