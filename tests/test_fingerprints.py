from uffpsim import mol_fp_supplier
import pytest

def test_calculate_fp_morgan():
    # Test invalid fp_type
    with pytest.raises(ValueError):
        list(mol_fp_supplier("data/10mols.smi", "InvalidFPType"))

    # Test valid input with default parameters
    fp_list_radius_2 = list(mol_fp_supplier("tests/data/10mols.smi", "Morgan"))
    assert isinstance(fp_list_radius_2, list)
    assert len(fp_list_radius_2) == 10
    assert all(isinstance(c, tuple) for c in fp_list_radius_2)
    assert all(len(c) == 3 for c in fp_list_radius_2)
    assert all(isinstance(c[0], str) and isinstance(c[1], str) and isinstance(c[2], str) for c in fp_list_radius_2)
    assert all(not c[1].isalpha() for c in fp_list_radius_2)
    assert all(len(c[1]) == 2048 for c in fp_list_radius_2)
    assert all(c[1].isdigit() for c in fp_list_radius_2)
    assert all(set(list(c[1])) == set(['0', '1']) for c in fp_list_radius_2)

    # Test valid input with custom parameters
    fp_list_radius_3 = list(mol_fp_supplier("tests/data/10mols.smi", "Morgan", {"radius": 3}))
    assert all(c2[1] != c3[1] for c2, c3 in zip(fp_list_radius_2, fp_list_radius_3))

    # Test valid input with custom parameters
    fp_list_bits_1024 = list(mol_fp_supplier("tests/data/10mols.smi", "Morgan", {"fpSize": 1024}))
    assert all(len(c[1]) == 1024 for c in fp_list_bits_1024)


def test_calculate_fp_avalon():
    # Test invalid fp_type
    with pytest.raises(ValueError):
        list(mol_fp_supplier("tests/data/10mols.smi", "InvalidFPType"))

    # Test valid input with default parameters
    fp_list = list(mol_fp_supplier("tests/data/10mols.smi", "Avalon"))
    assert isinstance(fp_list, list)
    assert len(fp_list) == 10
    assert all(isinstance(c, tuple) for c in fp_list)
    assert all(len(c) == 3 for c in fp_list)
    assert all(isinstance(c[0], str) and isinstance(c[1], str)  and isinstance(c[2], str) for c in fp_list)

def test_calculate_fp_topological():
    # Test invalid fp_type
    with pytest.raises(ValueError):
        list(mol_fp_supplier("tests/data/10mols.smi", "InvalidFPType"))

    # Test valid input with default parameters
    fp_list = list(mol_fp_supplier("tests/data/10mols.smi", "TopologicalTorsion"))
    assert isinstance(fp_list, list)
    assert len(fp_list) == 10
    assert all(isinstance(c, tuple) for c in fp_list)
    assert all(len(c) == 3 for c in fp_list)
    assert all(isinstance(c[0], str) and isinstance(c[1], str) and isinstance(c[2], str) for c in fp_list)
    assert all(not c[1].isalpha() for c in fp_list)
    assert all(len(c[1]) == 2048 for c in fp_list)

    # Test valid input with custom parameters
    fp_list_bits_1024 = list(mol_fp_supplier("tests/data/10mols.smi", "TopologicalTorsion", {"fpSize": 1024}))
    assert all(len(c[1]) == 1024 for c in fp_list_bits_1024)

def test_calculate_fp_rdkit():
    # Test invalid fp_type
    with pytest.raises(ValueError):
        list(mol_fp_supplier("tests/data/10mols.smi", "InvalidFPType"))

    # Test valid input with default parameters
    fp_list = list(mol_fp_supplier("tests/data/10mols.smi", "RDKit"))
    assert isinstance(fp_list, list)
    assert len(fp_list) == 10
    assert all(isinstance(c, tuple) for c in fp_list)
    assert all(len(c) == 3 for c in fp_list)
    assert all(isinstance(c[0], str) and isinstance(c[1], str) and isinstance(c[2], str) for c in fp_list)
    assert all(not c[1].isalpha() for c in fp_list)
    assert all(len(c[1]) == 2048 for c in fp_list)

    # Test valid input with custom parameters maxPath
    fp_list_max_path = list(mol_fp_supplier("tests/data/10mols.smi", "RDKit", {"maxPath": 6}))
    assert all(c2[1] != c3[1] for c2, c3 in zip(fp_list, fp_list_max_path))

    # Test valid input with custom parameters maxPath
    fp_list_bits_1024 = list(mol_fp_supplier("tests/data/10mols.smi", "RDKit", {"fpSize": 1024}))
    assert all(len(c[1]) == 1024 for c in fp_list_bits_1024)

def test_calculate_fp_atom_pair():
    # Test invalid fp_type
    with pytest.raises(ValueError):
        list(mol_fp_supplier("tests/data/10mols.smi", "InvalidFPType"))

    # Test valid input with default parameters
    fp_list = list(mol_fp_supplier("tests/data/10mols.smi", "AtomPair"))
    assert isinstance(fp_list, list)
    assert len(fp_list) == 10
    assert all(isinstance(c, tuple) for c in fp_list)
    assert all(len(c) == 3 for c in fp_list)
    assert all(isinstance(c[0], str) and isinstance(c[1], str)  and isinstance(c[2], str) for c in fp_list)
    assert all(not c[1].isalpha() for c in fp_list)
    assert all(len(c[1]) == 2048 for c in fp_list)

    # Test valid input with custom parameters maxDistance
    fp_list_max_dist = list(mol_fp_supplier("tests/data/10mols.smi", "AtomPair", {"maxDistance": 5}))
    assert all(c2[1] != c3[1] for c2, c3 in zip(fp_list, fp_list_max_dist))

    # Test valid input with custom parameters maxDistance
    fp_list_bits_1024 = list(mol_fp_supplier("tests/data/10mols.smi", "AtomPair", {"fpSize": 1024}))
    assert all(len(c[1]) == 1024 for c in fp_list_bits_1024)
