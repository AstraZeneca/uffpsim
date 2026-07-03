import timeit
import sys

from FPSim2 import FPSim2Engine

h5file=sys.argv[1]

smiles = []
with open("zinc_query_smiles.txt") as f:
    for line in f:
        smiles.append(line.split()[0])

engine = FPSim2Engine(h5file)

results = None
def sim_function(s):
    global results
    results = engine.top_k(s, threshold=0.7, k=1, metric='tanimoto', n_workers=1)

def run():
    for i, s in enumerate(smiles):
        spent_times = timeit.repeat('sim_function(s)', repeat=10, number=1, globals={**globals(), **locals()})

        if (len(results) > 0):
            print(i, 1, results[0][0], results[0][1], min(spent_times)*1000)
        else:
            print(i, 0, 0, 0, min(spent_times)*1000)

        #if i > 1:
        #    break


run()
