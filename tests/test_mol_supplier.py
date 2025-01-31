from uffpsim.fp_supplier import (
    smi_mol_supplier,
    sdf_mol_supplier,
    it_mol_supplier,
    get_mol_supplier
)
from rdkit import Chem
import os
import pytest

TESTS_DIR = os.path.dirname(os.path.abspath(__file__))

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
def test_sdf_mol_supplier_functionality():
    mol_supplier = sdf_mol_supplier("tests/data/10mols.sdf", gen_ids=True)
    mol_ids, mols, smiles = zip(*list(mol_supplier))
    assert mol_ids == ("1", "2", "3", "4", "5", "6", "7", "8", "9", "10")
    assert mols
    assert all(isinstance(mol, Chem.rdchem.Mol) for mol in mols)

    mol_supplier = sdf_mol_supplier("tests/data/10mols.sdf", gen_ids=False, mol_id_prop="mol_id")
    mol_ids, mols, smiles = zip(*list(mol_supplier))
    assert mol_ids ==  ("ID1", "ID2", "ID3", "ID4", "ID5", "ID6", "ID7", "ID8", "ID9", "ID10")
    assert mols
    assert all(isinstance(mol, Chem.rdchem.Mol) for mol in mols)

    mol_supplier = sdf_mol_supplier("tests/data/10mols.sdf", gen_ids=True, mol_id_prop="mol_id")
    mol_ids, mols, smiles = zip(*list(mol_supplier))
    assert mol_ids == ("1", "2", "3", "4", "5", "6", "7", "8", "9", "10")
    assert mols
    assert all(isinstance(mol, Chem.rdchem.Mol) for mol in mols)

    mol_supplier = sdf_mol_supplier("tests/data/10mols.sdf", gen_ids=False, mol_id_prop="not_a_mol_id")
    mol_ids, mols, smiles = zip(*list(mol_supplier))
    assert mol_ids == ("1", "2", "3", "4", "5", "6", "7", "8", "9", "10")
    assert mols
    assert all(isinstance(mol, Chem.rdchem.Mol) for mol in mols)

    mol_supplier = sdf_mol_supplier("tests/data/10mols.sdf.gz", gen_ids=False, mol_id_prop="mol_id")
    mol_ids, mols, smiles = zip(*list(mol_supplier))
    assert mol_ids ==  ("ID1", "ID2", "ID3", "ID4", "ID5", "ID6", "ID7", "ID8", "ID9", "ID10")
    assert mols
    assert all(isinstance(mol, Chem.rdchem.Mol) for mol in mols)

    mol_supplier = sdf_mol_supplier("tests/data/10mols.sdf.gz", gen_ids=True, mol_id_prop="mol_id")
    mol_ids, mols, smiles = zip(*list(mol_supplier))
    assert mol_ids == ("1", "2", "3", "4", "5", "6", "7", "8", "9", "10")
    assert mols
    assert all(isinstance(mol, Chem.rdchem.Mol) for mol in mols)

    mol_supplier = sdf_mol_supplier("tests/data/10mols.sdf.gz", gen_ids=False, mol_id_prop="not_a_mol_id")
    mol_ids, mols, smiles = zip(*list(mol_supplier))
    assert mol_ids == ("1", "2", "3", "4", "5", "6", "7", "8", "9", "10")
    assert mols
    assert all(isinstance(mol, Chem.rdchem.Mol) for mol in mols)

def test_smi_mol_supplier_functionality():
    mol_supplier = smi_mol_supplier("tests/data/10mols.smi", gen_ids=True)
    mol_ids, mols, smiles = zip(*list(mol_supplier))
    assert mol_ids == ("1", "2", "3", "4", "5", "6", "7", "8", "9", "10")
    assert mols
    assert all(isinstance(mol, Chem.rdchem.Mol) for mol in mols)

    mol_supplier = smi_mol_supplier("tests/data/10mols-with-ids.smi", gen_ids=True)
    mol_ids, mols, smiles = zip(*list(mol_supplier))
    assert mol_ids == ("1", "2", "3", "4", "5", "6", "7", "8", "9", "10")
    assert mols
    assert all(isinstance(mol, Chem.rdchem.Mol) for mol in mols)

    mol_supplier = smi_mol_supplier("tests/data/10mols-with-ids.smi", gen_ids=False)
    mol_ids, mols, smiles = zip(*list(mol_supplier))
    assert mol_ids == ("ID1", "ID2", "ID3", "ID4", "ID5", "ID6", "ID7", "ID8", "ID9", "ID10")
    assert mols
    assert all(isinstance(mol, Chem.rdchem.Mol) for mol in mols)

def test_it_mol_supplier_functionality():
    mol_supplier = it_mol_supplier(smiles_list, gen_ids=True)
    mol_ids, mols, smiles = zip(*list(mol_supplier))
    assert mol_ids == ("1", "2", "3", "4", "5", "6", "7", "8", "9", "10")
    assert mols
    assert all(isinstance(mol, Chem.rdchem.Mol) for mol in mols)

    mol_supplier = it_mol_supplier(smiles_list, gen_ids=False)
    mol_ids, mols, smiles = zip(*list(mol_supplier))
    assert mol_ids == ("1", "2", "3", "4", "5", "6", "7", "8", "9", "10")
    assert mols
    assert all(isinstance(mol, Chem.rdchem.Mol) for mol in mols)

    mol_supplier = it_mol_supplier([(s, f"ID{i+1}") for i,s in enumerate(smiles_list)], gen_ids=False)
    mol_ids, mols, smiles = zip(*list(mol_supplier))
    assert mol_ids == ("ID1", "ID2", "ID3", "ID4", "ID5", "ID6", "ID7", "ID8", "ID9", "ID10")
    assert mols
    assert all(isinstance(mol, Chem.rdchem.Mol) for mol in mols)

    mol_supplier = it_mol_supplier([(s, f"ID{i+1}") for i,s in enumerate(smiles_list)], gen_ids=True)
    mol_ids, mols, smiles = zip(*list(mol_supplier))
    assert mol_ids == ("1", "2", "3", "4", "5", "6", "7", "8", "9", "10")
    assert mols
    assert all(isinstance(mol, Chem.rdchem.Mol) for mol in mols)

def test_get_mol_supplier():
    assert get_mol_supplier("aaa.sdf") == sdf_mol_supplier
    assert get_mol_supplier("aaa.sdf.gz") == sdf_mol_supplier
    assert get_mol_supplier("aaa.smi") == smi_mol_supplier
    assert get_mol_supplier(smiles_list) == it_mol_supplier
