import glob
import numpy as np
import sys

def get_data():
    mol_data = []
    for filename in glob.glob(f'./logs/*.log',):
        with open(filename) as f:
            f.readline()
            for line in f:
                yield line

def get_selected_data(score_range):
    mol_data = [[], [], [], [], [], [], [], [], [], []]
    for data in get_data():
        d = data.split()
        score = np.round(float(d[3]), 3)
        popcount = int(d[2])
        if score >= score_range[0] and score < score_range[1] and popcount <= 10:
            mol_data[0].append(data)
        elif score >= score_range[0] and score < score_range[1] and popcount <= 20:
            mol_data[1].append(data)
        elif score >= score_range[0] and score < score_range[1] and popcount <= 30:
            mol_data[2].append(data)
        elif score >= score_range[0] and score < score_range[1] and popcount <= 40:
            mol_data[3].append(data)
        elif score >= score_range[0] and score < score_range[1] and popcount <= 50:
            mol_data[4].append(data)
        elif score >= score_range[0] and score < score_range[1] and popcount <= 60:
            mol_data[5].append(data)
        elif score >= score_range[0] and score < score_range[1] and popcount <= 70:
            mol_data[6].append(data)
        elif score >= score_range[0] and score < score_range[1] and popcount <= 80:
            mol_data[7].append(data)
        elif score >= score_range[0] and score < score_range[1] and popcount <= 90:
            mol_data[8].append(data)
        elif score >= score_range[0] and score < score_range[1] and popcount > 90:
            mol_data[9].append(data)

    selected_data = []
    for i, data in enumerate(mol_data):
        print(i, len(data))
        if len(data) > 0:
            last = 500
            if len(data) < 500:
                last = len(data) - 1
            idx = np.random.choice(len(data), last, replace=False)
            np_data = np.array(data)
            s_data = np_data[idx]
            selected_data += list(s_data)

    return selected_data


with open('selected_data.txt', 'w') as f:
    for srange in [[0.5, 0.6], [0.6, 0.7], [0.7, 0.8], [0.8, 0.9], [0.9, 0.99], [1.0, 1.1]]:
        slected_data = get_selected_data(srange)
        for data in slected_data:
            f.write(data)
