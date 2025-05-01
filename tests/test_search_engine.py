import os
from uffpsim import create_database, UFFPSimSearchEngine

def test_search_sequential():
    db_file = 'test.h5'
    info = {"createdAt": "25/01/2025", "version": 2025}
    fp_params = {"radius": 2, "fpSize": 2048}
    create_database("tests/data/10mols.sdf", db_file, "Morgan", fp_params, gen_ids=False, mol_id_prop="mol_id",
                     info=info, inner_clustering_threshold=0.1)
    
    search_engine = UFFPSimSearchEngine("test.h5")

    smiles_list = [
        "Cc1cc(-n2ncc(=O)[nH]c2=O)ccc1C(=O)c1ccccc1Cl",
        "Cc1cc(-n2ncc(=O)[nH]c2=O)ccc1C(=O)c1ccc(C#N)cc1",
        "Cc1cc(-n2ncc(=O)[nH]c2=O)cc(C)c1C(O)c1ccc(Cl)cc1",
        "Cc1ccc(C(=O)c2ccc(-n3ncc(=O)[nH]c3=O)cc2)cc1",
        "Cc1cc(-n2ncc(=O)[nH]c2=O)ccc1C(=O)c1ccc(Cl)cc1",
        "Cc1cc(-n2ncc(=O)[nH]c2=O)ccc1C(=O)c1ccccc1",
        "Cc1cc(Br)ccc1C(=O)c1ccc(-n2ncc(=O)[nH]c2=O)cc1Cl",
        "O=C(c1ccc(Cl)cc1Cl)c1ccc(-n2ncc(=O)[nH]c2=O)cc1Cl",
        "CS(=O)(=O)c1ccc(C(=O)c2ccc(-n3ncc(=O)[nH]c3=O)cc2Cl)cc1",
        "c1cc2cc(c1)-c1cccc(c1)C[n+]1ccc(c3ccccc31)NCCCCCCCCCCNc1cc[n+](c3ccccc13)C2",
    ]

    for idx, smiles in enumerate(smiles_list):
        r = search_engine.search(smiles, 0.8, 1)
        assert r is not None
        assert len(r) == 1
        assert r[0][0] == f"ID{idx+1}"

    # Remove the database file
    if os.path.isfile("test.h5"):
        os.remove("test.h5")

def test_search_bad_smiles_sequential():
    db_file = 'test.h5'
    info = {"createdAt": "25/01/2025", "version": 2025}
    fp_params = {"radius": 2, "fpSize": 2048}
    create_database("tests/data/10mols.sdf", db_file, "Morgan", fp_params, gen_ids=False, mol_id_prop="mol_id",
                     info=info, inner_clustering_threshold=0.1)
    
    search_engine = UFFPSimSearchEngine("test.h5")

    smiles_list = [
        "Cc1cc(-n2ncc(=O)[nH]c2=O)ccc1C(=O)c1ccccc1Cl",
        "Cc1cc(-n2ncc(=O)[nH]c2=O)ccc1C(=O)c1ccc(C#N)cc1",
        "Cc1cc(-n2ncc(=O)[nH]c2=O)cc(C)c1C(O)c1ccc(Cl)cc1",
        "Cc1ccc(C(=O)c2ccc(-n3ncc)[nH]c3=O)cc2)cc1",   # bad smiles, result should be None
        "Cc1cc(-n2ncc(=O)[nH]c2=O)ccc1C(=O)c1ccc(Cl)cc1",
        "Cc1cc(-n2ncc(=O)[nH]c2=O)ccc1C(=O)c1ccccc1",
        "Cc1cc(Br)ccc1C(=O)c1ccc(-n2ncc(=O)[nH]c2=O)cc1Cl",
        "O=C(c1ccc(Cl)cc1Cl)c1ccc=O)[nH]c2=O)cc1Cl",   # bad smiles, result should be None
        "CS(=O)(=O)c1ccc(C(=O)c2ccc(-n3ncc(=O)[nH]c3=O)cc2Cl)cc1",
        "c1cc2cc(c1)-c1cccc(c1)C[n+]1ccc(c3ccccc31)NCCCCCCCCCCNc1cc[n+](c3ccccc13)C2",
    ]

    for idx, smiles in enumerate(smiles_list):
        r = search_engine.search(smiles, 0.8, 1)
        if idx == 3 or idx == 7:
            assert r is None
        else:
            assert r is not None
            assert len(r) == 1
            assert r[0][0] == f"ID{idx+1}"

    # Remove the database file
    if os.path.isfile("test.h5"):
        os.remove("test.h5")

def test_search_batch():
    db_file = 'test.h5'
    info = {"createdAt": "25/01/2025", "version": 2025}
    fp_params = {"radius": 2, "fpSize": 2048}
    create_database("tests/data/10mols.sdf", db_file, "Morgan", fp_params, gen_ids=False, mol_id_prop="mol_id",
                     info=info, inner_clustering_threshold=0.1)
    
    search_engine = UFFPSimSearchEngine("test.h5")

    smiles_list = [
        "Cc1cc(-n2ncc(=O)[nH]c2=O)ccc1C(=O)c1ccccc1Cl",
        "Cc1cc(-n2ncc(=O)[nH]c2=O)ccc1C(=O)c1ccc(C#N)cc1",
        "Cc1cc(-n2ncc(=O)[nH]c2=O)cc(C)c1C(O)c1ccc(Cl)cc1",
        "Cc1ccc(C(=O)c2ccc(-n3ncc(=O)[nH]c3=O)cc2)cc1",
        "Cc1cc(-n2ncc(=O)[nH]c2=O)ccc1C(=O)c1ccc(Cl)cc1",
        "Cc1cc(-n2ncc(=O)[nH]c2=O)ccc1C(=O)c1ccccc1",
        "Cc1cc(Br)ccc1C(=O)c1ccc(-n2ncc(=O)[nH]c2=O)cc1Cl",
        "O=C(c1ccc(Cl)cc1Cl)c1ccc(-n2ncc(=O)[nH]c2=O)cc1Cl",
        "CS(=O)(=O)c1ccc(C(=O)c2ccc(-n3ncc(=O)[nH]c3=O)cc2Cl)cc1",
        "c1cc2cc(c1)-c1cccc(c1)C[n+]1ccc(c3ccccc31)NCCCCCCCCCCNc1cc[n+](c3ccccc13)C2",
    ]

    result = search_engine.batch_search(smiles_list, 0.8, 1)
    assert result is not None
    assert len(result) == 10
    for idx, r in enumerate(result):
        assert len(r) == 1
        assert r[0][0] == f"ID{idx+1}"

    # Remove the database file
    if os.path.isfile("test.h5"):
        os.remove("test.h5")

def test_search_bad_smiles_batch():
    db_file = 'test.h5'
    info = {"createdAt": "25/01/2025", "version": 2025}
    fp_params = {"radius": 2, "fpSize": 2048}
    create_database("tests/data/10mols.sdf", db_file, "Morgan", fp_params, gen_ids=False, mol_id_prop="mol_id",
                     info=info, inner_clustering_threshold=0.1)
    
    search_engine = UFFPSimSearchEngine("test.h5")

    smiles_list = [
        "Cc1cc(-n2ncc(=O)[nH]c2=O)ccc1C(=O)c1ccccc1Cl",
        "Cc1cc(-n2ncc(=O)[nH]c2=O)ccc1C(=O)c1ccc(C#N)cc1",
        "Cc1cc(-n2ncc(=O)[nH]c2=O)cc(C)c1C(O)c1ccc(Cl)cc1",
        "Cc1ccc(C(=O)c2ccc(-n3ncc)[nH]c3=O)cc2)cc1",   # bad smiles, result should be None
        "Cc1cc(-n2ncc(=O)[nH]c2=O)ccc1C(=O)c1ccc(Cl)cc1",
        "Cc1cc(-n2ncc(=O)[nH]c2=O)ccc1C(=O)c1ccccc1",
        "Cc1cc(Br)ccc1C(=O)c1ccc(-n2ncc(=O)[nH]c2=O)cc1Cl",
        "O=C(c1ccc(Cl)cc1Cl)c1ccc=O)[nH]c2=O)cc1Cl",   # bad smiles, result should be None
        "CS(=O)(=O)c1ccc(C(=O)c2ccc(-n3ncc(=O)[nH]c3=O)cc2Cl)cc1",
        "c1cc2cc(c1)-c1cccc(c1)C[n+]1ccc(c3ccccc31)NCCCCCCCCCCNc1cc[n+](c3ccccc13)C2",
    ]

    result = search_engine.batch_search(smiles_list, 0.8, 1)
    assert result is not None
    assert len(result) == 10
    for idx, r in enumerate(result):
        if idx == 3 or idx == 7:
            assert r is None
        else:
            assert len(r) == 1
            assert r[0][0] == f"ID{idx+1}"

    # Remove the database file
    if os.path.isfile("test.h5"):
        os.remove("test.h5")

def test_search_smiles_for_mol_id():
    db_file = 'test.h5'
    info = {"createdAt": "25/01/2025", "version": 2025}
    fp_params = {"radius": 2, "fpSize": 2048}
    create_database("tests/data/10mols.sdf", db_file, "Morgan", fp_params, gen_ids=False, mol_id_prop="mol_id",
                     info=info, inner_clustering_threshold=0.1)
    
    search_engine = UFFPSimSearchEngine("test.h5")

    smiles_list = [
        "Cc1cc(-n2ncc(=O)[nH]c2=O)ccc1C(=O)c1ccccc1Cl",
        "Cc1cc(-n2ncc(=O)[nH]c2=O)ccc1C(=O)c1ccc(C#N)cc1",
        "Cc1cc(-n2ncc(=O)[nH]c2=O)cc(C)c1C(O)c1ccc(Cl)cc1",
        "Cc1ccc(C(=O)c2ccc(-n3ncc(=O)[nH]c3=O)cc2)cc1",
        "Cc1cc(-n2ncc(=O)[nH]c2=O)ccc1C(=O)c1ccc(Cl)cc1",
        "Cc1cc(-n2ncc(=O)[nH]c2=O)ccc1C(=O)c1ccccc1",
        "Cc1cc(Br)ccc1C(=O)c1ccc(-n2ncc(=O)[nH]c2=O)cc1Cl",
        "O=C(c1ccc(Cl)cc1Cl)c1ccc(-n2ncc(=O)[nH]c2=O)cc1Cl",
        "CS(=O)(=O)c1ccc(C(=O)c2ccc(-n3ncc(=O)[nH]c3=O)cc2Cl)cc1",
        "c1cc2cc(c1)-c1cccc(c1)C[n+]1ccc(c3ccccc31)NCCCCCCCCCCNc1cc[n+](c3ccccc13)C2",
    ]

    result = search_engine.batch_search(smiles_list, 0.8, 1)
    assert result is not None
    assert len(result) == 10

    for idx, r in enumerate(result):
        assert len(r) == 1
        assert r[0][0] == f"ID{idx+1}"
        mol_id = f"ID{idx+1}"
        smiles = search_engine.get_smiles_for_id(mol_id)
        assert smiles == smiles_list[idx]