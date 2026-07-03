"""Download ChEMBL 33 SDF and create the UFFPSim database used by 01-screen.py."""
import subprocess
from pathlib import Path

from uffpsim import create_database


CHEMBL_URL = "https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases/chembl_33/chembl_33.sdf.gz"
BASE_DIR = Path(__file__).parent / "../.."
SDF_FILE = BASE_DIR / "chembl_33.sdf.gz"
DB_FILE = BASE_DIR / "chembl_2048b.h5"


def download_chembl():
    if SDF_FILE.exists():
        print(f"{SDF_FILE} already exists, skipping download.")
        return

    print(f"Downloading {CHEMBL_URL} ...")
    result = subprocess.run(
        ["curl", "--fail", "--location", "--progress-bar", "-o", str(SDF_FILE), CHEMBL_URL],
        check=True,
    )
    print(f"Downloaded to {SDF_FILE}")


def build_database():
    if DB_FILE.exists():
        print(f"{DB_FILE} already exists, skipping database creation.")
        return

    print(f"Creating database {DB_FILE} ...")
    create_database(
        input_file=str(SDF_FILE),
        db_file=str(DB_FILE),
        fp_type="Morgan",
        fp_params={"fpSize": 2048, "radius": 2},
        gen_ids=True,
        mol_id_prop=None,
        mol_id_max_chars=15,
        inner_clustering_threshold=0.15,
    )
    print(f"Database created: {DB_FILE}")


if __name__ == "__main__":
    download_chembl()
    build_database()
