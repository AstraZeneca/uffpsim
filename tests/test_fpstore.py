import os
import json
from uffpsim import create_database, update_database, uffpsimLib, redo_inner_clustering

def test_create_database():
    # Test valid input file
    info = {"createdAt": "25/01/2025", "version": 2025}
    fp_params = {"radius": 2, "fpSize": 2048}
    create_database("tests/data/10mols.sdf", "test.h5", "Morgan", fp_params, gen_ids=True, mol_id_prop=None,
                     info=info, inner_clustering_threshold=0.1)

    # Test if database file is created
    assert os.path.isfile("test.h5")

    fp_store = uffpsimLib.FingerprintStore("test.h5")


    # check fingerprint parameters stored in file
    fp_params_from_file = json.loads(fp_store.fp_params_json)
    assert fp_params_from_file["fp_type"] == "Morgan"
    assert fp_params_from_file["fp_params"]["radius"] == 2
    assert fp_params_from_file["fp_params"]["fpSize"] == 2048

    # check attributes of FingerprintStore read from file
    assert fp_store.mol_id_max_chars == 15
    assert int(fp_store.inner_clustering_threshold * 10) == 1
    assert fp_store.fp_bits_size == 2048
    assert fp_store.num_mols == 10

    #check info
    assert json.loads(fp_store.info) == info

    # check smiles from id
    assert fp_store.get_smiles_for_id("1") == 'Cc1cc(-n2ncc(=O)[nH]c2=O)ccc1C(=O)c1ccccc1Cl'
    assert fp_store.get_smiles_for_id("5") == 'Cc1cc(-n2ncc(=O)[nH]c2=O)ccc1C(=O)c1ccc(Cl)cc1'
    assert fp_store.get_smiles_for_id("10") == 'c1cc2cc(c1)-c1cccc(c1)C[n+]1ccc(c3ccccc31)NCCCCCCCCCCNc1cc[n+](c3ccccc13)C2'

    fp_store.close()

    # Remove the database file
    if os.path.isfile("test.h5"):
        os.remove("test.h5")

def test_redo_inner_clustering():
    db_file = 'test1.h5'
    info = {"createdAt": "25/01/2025", "version": 2025}
    fp_params = {"radius": 2, "fpSize": 2048}
    create_database("tests/data/10mols.sdf", db_file, "Morgan", fp_params, gen_ids=True, mol_id_prop=None,
                     info=info, inner_clustering_threshold=0.1)

    # Test if database file is created
    assert os.path.isfile(db_file)

    fp_store = uffpsimLib.FingerprintStore("test1.h5")
    assert int(fp_store.inner_clustering_threshold*10) == 1

    fp_store.close()

    # Test redoing inner clustering
    redo_inner_clustering(db_file, 0.85)
    fp_store = uffpsimLib.FingerprintStore(db_file)
    assert int(fp_store.inner_clustering_threshold*100) == 85
    fp_store.close()

    # Remove the database file
    if os.path.isfile(db_file):
        os.remove(db_file)

def test_mol_id_clustering():
    db_file = 'test2.h5'
    info = {"createdAt": "25/01/2025", "version": 2025}
    fp_params = {"radius": 2, "fpSize": 2048}
    create_database("tests/data/10mols.sdf", db_file, "Morgan", fp_params, gen_ids=False, mol_id_prop="mol_id",
                     info=info, inner_clustering_threshold=0.1)
    
    # Test if database file is created
    assert os.path.isfile(db_file)

    fp_store = uffpsimLib.FingerprintStore(db_file)

    # check smiles from id
    assert fp_store.get_smiles_for_id("ID1") == 'Cc1cc(-n2ncc(=O)[nH]c2=O)ccc1C(=O)c1ccccc1Cl'
    assert fp_store.get_smiles_for_id("ID5") == 'Cc1cc(-n2ncc(=O)[nH]c2=O)ccc1C(=O)c1ccc(Cl)cc1'
    assert fp_store.get_smiles_for_id("ID10") == 'c1cc2cc(c1)-c1cccc(c1)C[n+]1ccc(c3ccccc31)NCCCCCCCCCCNc1cc[n+](c3ccccc13)C2'

    fp_store.close()

    # Remove the database file
    if os.path.isfile(db_file):
        os.remove(db_file)

def test_update_fpstore():
    db_file = 'test3.h5'
    info = {"createdAt": "25/01/2025", "version": 2025}
    fp_params = {"radius": 2, "fpSize": 2048}
    create_database("tests/data/10mols.sdf", db_file, "Morgan", fp_params, gen_ids=False, mol_id_prop="mol_id",
                     info=info, inner_clustering_threshold=0.1)
    
    # Test if database file is created
    assert os.path.isfile(db_file)

    fp_store = uffpsimLib.FingerprintStore(db_file)

    # check smiles from id
    assert fp_store.get_smiles_for_id("ID1") == 'Cc1cc(-n2ncc(=O)[nH]c2=O)ccc1C(=O)c1ccccc1Cl'
    assert fp_store.get_smiles_for_id("ID3") == 'Cc1cc(-n2ncc(=O)[nH]c2=O)cc(C)c1C(O)c1ccc(Cl)cc1'

    fp_store.close()

    smiles_list = [
        ("Cc1cc(-n2ncc(=O)[nH]c2=O)ccc1C(=O)c1ccncc1Cl", "ID1"),
        ("Cc1cc(-n2ncc(=O)[nH]c2=O)cc(C)c1C(O)c1cnc(Cl)cc1", "ID3")
    ]

    update_database(smiles_list, db_file=db_file)

    fp_store = uffpsimLib.FingerprintStore(db_file)

    # check updated smiles from id
    assert fp_store.get_smiles_for_id("ID1") == 'Cc1cc(-n2ncc(=O)[nH]c2=O)ccc1C(=O)c1ccncc1Cl'
    assert fp_store.get_smiles_for_id("ID3") == 'Cc1cc(-n2ncc(=O)[nH]c2=O)cc(C)c1C(O)c1cnc(Cl)cc1'

    fp_store.close()
    
    # Remove the database file
    if os.path.isfile(db_file):
        os.remove(db_file)