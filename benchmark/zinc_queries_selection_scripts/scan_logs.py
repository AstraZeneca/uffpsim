import glob

import numpy as np


def get_data():
    for filename in glob.glob("./logs/*.log"):
        with open(filename) as f:
            f.readline()
            for line in f:
                yield line


def get_selected_data(score_range):
    mol_data = [[], [], [], [], [], [], [], [], [], []]

    for data in get_data():
        fields = data.split()
        popcount = int(fields[2])
        score = np.round(float(fields[3]), 3)

        if not (score_range[0] <= score < score_range[1]):
            continue

        if popcount <= 10:
            mol_data[0].append(data)
        elif popcount <= 20:
            mol_data[1].append(data)
        elif popcount <= 30:
            mol_data[2].append(data)
        elif popcount <= 40:
            mol_data[3].append(data)
        elif popcount <= 50:
            mol_data[4].append(data)
        elif popcount <= 60:
            mol_data[5].append(data)
        elif popcount <= 70:
            mol_data[6].append(data)
        elif popcount <= 80:
            mol_data[7].append(data)
        elif popcount <= 90:
            mol_data[8].append(data)
        else:
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
            selected_data += list(np_data[idx])

    return selected_data


with open("selected_data.txt", "w") as f:
    for score_range in [
        [0.5, 0.6],
        [0.6, 0.7],
        [0.7, 0.8],
        [0.8, 0.9],
        [0.9, 0.99],
        [1.0, 1.1],
    ]:
        selected_data = get_selected_data(score_range)

        for data in selected_data:
            f.write(data)
