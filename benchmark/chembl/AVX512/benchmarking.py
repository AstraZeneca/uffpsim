import json
import csv
import timeit
import pickle

from uffpsim import UFFPSimSearchEngine

h5file='chembl_2048b.h5'

smiles = []
with open("zinc_query_smiles.txt") as f:
    for line in f:
        smiles.append(line.split()[0])

engine = UFFPSimSearchEngine(h5file)

results = None
def sim_function(s):
    global results
    results = engine.search(s, 0.7, 1)

def run():
    for i, s in enumerate(smiles):
        spent_times = timeit.repeat('sim_function(s)', repeat=20, number=1, globals={**globals(), **locals()})

        if (len(results) > 0):
            print(i, 1, results[0][0], results[0][1], min(spent_times)*1000)
        else:
            print(i, 0, 0, 0, min(spent_times)*1000)

        #if i > 1:
        #    break


run()
