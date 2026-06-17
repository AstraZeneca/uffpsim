
from FPSim2.io import create_db_file

for size in [256, 512, 1024, 2048]:
    create_db_file(
    mols_source='../chembl_35.smi',
    filename=f'chembl_35_b{size}.h5',
    mol_format=None, # set to None
    fp_type='Morgan',
    fp_params={'radius': 2, 'fpSize': size},
    mol_id_prop=None
    )
