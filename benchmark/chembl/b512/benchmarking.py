import json
import csv
import timeit
import pickle
import sys

from uffpsim import UFFPSimSearchEngine

h5file = sys.argv[1]

smiles = []
with open("../../zinc_query_smiles.txt") as f:
    for line in f:
        smiles.append(line.split()[0])

engine = UFFPSimSearchEngine(h5file)

results = None
def sim_function(s):
    global results
    results = engine.search(s, 0.7, 1)

def run():
    for i, s in enumerate(smiles):
        spent_times = timeit.repeat('sim_function(s)', repeat=10, number=1, globals={**globals(), **locals()})
        hits = results[0]
        ops = results[1]
        if (len(hits) > 0):
            print(i, 1, hits[0][0], hits[0][1], min(spent_times)*1000, ops)
        else:
            print(i, 0, 0, 0, min(spent_times)*1000, ops)

        #if i > 1:
        #    break


run()
