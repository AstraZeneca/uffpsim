from uffpsim import redo_inner_clustering
import shutil
import os

cutoffs =  {"10": 0.1, "11": 0.11, "12": 0.12,"13": 0.13, "14": 0.14,
            "16": 0.16, "17": 0.17, "18": 0.18, "19": 0.19, "20": 0.2,
            "21": 0.21, "22": 0.22, "23": 0.23, "24": 0.24, "25": 0.25,
            "26": 0.26, "27": 0.27, "28": 0.28, "29": 0.29, "30": 0.30,
            "31": 0.31, "32": 0.32, "33": 0.33, "34": 0.34, "35": 0.35}

if not os.path.exists("chembl35_20.h5"):
    raise FileNotFoundError("chembl35_20.h5 not found. Please run create_db.py first to generate the initial database.")

for key, value in cutoffs.items():
    shutil.copyfile("chembl35_20.h5", f"chembl35_{key}.h5")
    redo_inner_clustering(f"chembl35_{key}.h5", value, cluster_mode= "memory", cluster_parallel=False)
