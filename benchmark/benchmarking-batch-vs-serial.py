import json
import csv
import timeit
import numpy as np
import pickle

from uffpsim import UFFPSimSearchEngine

h5file='chembl_2048b.h5'

smiles = []
with open("zinc_query_smiles.txt") as f:
    for line in f:
        smiles.append(line.split()[0])

engine = UFFPSimSearchEngine(h5file)

def sim_function_batch(smlist):
    engine.batch_search(smlist, 0.7, 1)

results = None
def sim_function_serial(s):
    engine.search(s, 0.7, 1)

def run():
    for iter, sample_size in [(1000, 100), (200, 500), (100, 1000), (50, 2000), (25, 4000)]:
        for i in range(iter):
            smlist = list(np.random.choice(smiles, sample_size, replace=False))
            batch_spent_times = timeit.repeat('sim_function_batch(smlist)', repeat=5, number=1, globals={**globals(), **locals()})
            batch_spent_times = min(batch_spent_times)*1000
            serial_spent_times = 0
            for s in smlist:
                tmp_spent_times = timeit.repeat('sim_function_serial(s)', repeat=5, number=1, globals={**globals(), **locals()})
                serial_spent_times += min(tmp_spent_times)*1000

            print(i, sample_size, batch_spent_times, serial_spent_times, serial_spent_times/batch_spent_times)

        #if i > 1:
        #    break


run()
