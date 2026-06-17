from uffpsim import create_database_parallel

info = {"downloaded": "06/03/2025", "version": "35"}

for size in [256, 512, 1024, 2048]:
    fp_params = {"radius": 2, "fpSize": size}
    
    create_database_parallel("chembl_35.sdf", f"b{size}/chembl35_20.h5", "Morgan", workers= 8, fp_params=fp_params, 
                    gen_ids=False, mol_id_prop="chembl_id", mol_id_max_chars = 15, info=info,
                    inner_clustering_threshold = 0.20, cluster_mode = "memory", cluster_parallel = True)
