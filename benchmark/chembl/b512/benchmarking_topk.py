import json
import csv
import timeit
import pickle
import sys

from uffpsim import UFFPSimSearchEngine

h5file = sys.argv[1]
cutoff = float(sys.argv[2])
topk = int(sys.argv[3])

smiles = []
with open("../../zinc_query_smiles.txt") as f:
    for line in f:
        smiles.append(line.split()[0])

engine = UFFPSimSearchEngine(h5file)

results = None
def sim_function(s):
    global results
    results = engine.search(s, cutoff, topk)

def run():
    for i, s in enumerate(smiles):
        spent_times = timeit.repeat('sim_function(s)', repeat=10, number=1, globals={**globals(), **locals()})
        hits = results[0]
        ops = results[1]
        if (len(hits) > 0):
            print(i, len(hits), hits[0][0], hits[0][1], hits[-1][0], hits[-1][1], min(spent_times)*1000, ops)
        else:
            print(i, 0, 0, 0, 0, 0, min(spent_times)*1000, ops)

        #if i > 1:
        #    break


run()
