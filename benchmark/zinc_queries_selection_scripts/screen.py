import glob
import sys
from uffpsim import UFFPSimSearchEngine

def get_smiles_batch():
    count = 0
    mol_data = []
    for filename in glob.glob(f'../data/{sys.argv[1]}/*/*.smi',):
        with open(filename) as f:
            f.readline()
            for line in f:
                d = line.split()
                mol_data.append(d)

                if len(mol_data) == 1000:
                    yield mol_data
                    mol_data = []

    yield mol_data




search_engine = UFFPSimSearchEngine("../../chembl_2048b.h5")

for mol_data in get_smiles_batch():
    result = search_engine.batch_search([d[0] for d in mol_data], 0.6, 1)
    for d, r in zip(mol_data, result):
        if len(r) > 0:
            smiles = d[0]
            zincId = d[1]
            cfp = search_engine.getCompactFingerPrintArray(smiles)
            print(smiles, zincId, cfp[-1], r[0][1])