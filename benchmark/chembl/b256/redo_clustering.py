from uffpsim import redo_inner_clustering
import shutil
import os

cutoffs =  {"10": 0.10, "15": 0.15, "10": 0.10, "20": 0.2, "22": 0.22, "24": 0.24, "26": 0.26, "28": 0.28, "30": 0.30,
            "32": 0.32, "34": 0.34, "36": 0.36, "38": 0.38, "40": 0.40, "42": 0.42, "44": 0.44}

if not os.path.exists("chembl35_20.h5"):
    raise FileNotFoundError("chembl35_20.h5 not found. Please run create_db.py first to generate the initial database.")

for key, value in cutoffs.items():
    shutil.copyfile(f"chembl35_20.h5", f"chembl35_{key}.h5")
    redo_inner_clustering(f"chembl35_{key}.h5", value, cluster_mode= "memory", cluster_parallel=True)
