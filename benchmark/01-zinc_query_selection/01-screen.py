from pathlib import Path
import sys

from uffpsim import UFFPSimSearchEngine


BATCH_SIZE = 1000
SEARCH_THRESHOLD = 0.5
TOP_K = 1


def get_smiles_batches(dataset: str, batch_size: int = BATCH_SIZE):
    batch = []

    for smi_file in Path("../data", dataset).glob("*/*.smi"):
        with smi_file.open() as f:
            next(f)  # Skip header

            for line in f:
                batch.append(line.split())

                if len(batch) >= batch_size:
                    yield batch
                    batch = []

    if batch:
        yield batch


def main():
    if len(sys.argv) != 2:
        print(f"Usage: {sys.argv[0]} <dataset>")
        sys.exit(1)

    search_engine = UFFPSimSearchEngine("../../chembl_2048b.h5")

    for batch in get_smiles_batches(sys.argv[1]):
        smiles_list = [smiles for smiles, *_ in batch]
        results = search_engine.batch_search(
            smiles_list,
            SEARCH_THRESHOLD,
            TOP_K,
        )

        for (smiles, zinc_id, *_), matches in zip(batch, results):
            if not matches:
                continue

            cfp = search_engine.getCompactFingerPrintArray(smiles)
            print(smiles, zinc_id, cfp[-1], matches[0][1])


if __name__ == "__main__":
    main()
