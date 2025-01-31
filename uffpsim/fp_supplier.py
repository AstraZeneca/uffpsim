"""
Fingerprint suppliers in serial. 
Taken from FPSim2 source code.

"""
from typing import Any, Callable, Iterable as IterableType, Dict, List, Tuple, Union
from collections.abc import Iterable
from rdkit import rdBase
from rdkit import Chem

from .fingerprints import FPCalculator, load_molecule

def mol_fp_supplier(filename: str, fp_type: str, fp_params: Dict[str, Any] = None, gen_ids: bool = True, mol_id_prop=None, mol_id_length: int = 15) -> List[Tuple[str, str, str]]:
    """Calculates the fingerprints for the molecules in the given file.
    The fingerprints are returned as a list of tuples (mol_id, fingerprint).
    The mol_id can be used to match the fingerprints with the original molecules.
    The fingerprints are returned as bit strings.
    The fingerprints are generated using the given fingerprint type and parameters.
    The default parameters are used if no parameters are given.

    Parameters
    ----------
    filename : str
        The path to the file containing the molecules.
    fp_type : str
        The type of fingerprint to calculate.
    fp_params : Dict[str, Any], optional
        A dictionary of parameters to use for the fingerprint calculation.
        The default is None.
    gen_ids : bool, optional
        Whether to automatically generate IDs for the molecules. If False, the
        IDs must be provided in the input file. The default is True.
    mol_id_prop : Union[str, None], optional
        The name of the property in the input file that contains the molecule
        IDs. If None, the IDs will be generated automatically. The default is None.

    Returns
    -------
    List[Tuple[str, str]]
        A list of tuples containing the molecule IDs and their fingerprints. The
        fingerprints are represented as bit strings.

    Raises
    ------
    ValueError
        If the fingerprint type is not supported.
    ValueError
        If an invalid parameter is provided for the fingerprint type.
    ValueError
        If the input file cannot be read.

    Examples
    --------
    >>> morgan_fingerprints = list(calculate_fp("tests/data/10mols.smi", "Morgan"))
    >>> # from list of smiles:
    >>> fingerprints = list(calculate_fp(["Cc1cc(-n2ncc(=O)[nH]c2=O)ccc1C(=O)c1ccccc1Cl", "O=C(c1ccc(Cl)cc1Cl)c1ccc(-n2ncc(=O)[nH]c2=O)cc1Cl"], "Morgan"))
    >>> # show default parameters for Morgan fingerprints:
    >>> print(uffpsim.fingerprints.FP_TYPE_DEFAULT_PARAMETERS["Morgan"])
    >>> # show all available fingerprint types:
    >>> print(uffpsim.fingerprints.FP_TYPES)
    """
    blocker = rdBase.BlockLogs()
    if (mol_id_length+1) % 8 != 0:
        raise ValueError(f"mol_id_length {mol_id_length} must be less than one multiple of 8. e.g. 7, 15, 23, 31 etc.")
    
    fp_calculator = FPCalculator(fp_type, fp_params)

    supplier = get_mol_supplier(filename)
    if not supplier:
        raise ValueError("Invalid input! Cannot find a mol supplier for the given input!")
    
    for mol_id, rdmol, smiles in supplier(filename, gen_ids=gen_ids, mol_id_prop=mol_id_prop):
        if len(mol_id) >= mol_id_length:
            raise ValueError(f"mol_id must be less than {mol_id_length} characters long! Found mol_id: {mol_id}!")
        
        if rdmol is not None:
            yield (mol_id, fp_calculator(rdmol), smiles)

def it_mol_supplier(iterable: IterableType, gen_ids: bool=True, **kwargs) -> IterableType[Tuple[int, Chem.Mol, str]]:
    """Generator function that reads from iterables.

    Parameters
    ----------
    iterable : iterable
         Python iterable storing molecules.

    gen_ids: bool
        if True, automatically generate ids otherwise use from input file.

    Yields
    -------
    tuple
        int id and rdkit mol.
    """
    for new_mol_id, mol in enumerate(iterable, 1):
        show_warning = True
        if isinstance(mol, str):
            mol_string = mol
            mol_id = new_mol_id
        else:
            if gen_ids:
                mol_string = mol[0]
                mol_id = new_mol_id
            else:
                if len(mol) == 1:
                    mol_string = mol
                    mol_id = new_mol_id
                    if show_warning:
                        show_warning = False
                        raise  Warning("No mol_id found in the iterable!"
                                        "Consider setting gen_ids=True when running to automatically generate the mol_id from index.")
                        
                else:
                    mol_string = mol[0]
                    mol_id = mol[1]
            
        (rdmol, smiles) = load_molecule(mol_string)
        if rdmol:
            yield str(mol_id), rdmol, smiles
        else:
            continue


def smi_mol_supplier(filename: str, gen_ids: bool=True, **kwargs) -> IterableType[Tuple[str, Chem.Mol, str]]:
    """Generator function that reads from a .smi file.

    Parameters
    ----------
    filename : str
         .smi file name.

    gen_ids: bool
        if True, automatically generate ids otherwise use from input file.

    Yields
    -------
    tuple
        int id and rdkit mol.
    """
    with open(filename, "r") as f:
        show_warning = True
        for new_mol_id, mol in enumerate(f, 1):
            mol = mol.split()
            if gen_ids:
                smiles_or_inchi = mol[0]
                mol_id = new_mol_id
            else:
                if len(mol) == 1:
                    smiles_or_inchi = mol[0]
                    mol_id = new_mol_id
                    if show_warning:
                        show_warning = False
                        print("WARNING: No mol_id found in the SMILES file!"
                                       "Consider setting gen_ids=True when running to automatically generate the mol_id from index.")
                else:
                    smiles_or_inchi = mol[0]
                    mol_id = mol[1]
            
            smiles_or_inchi = smiles_or_inchi.strip()
            if smiles_or_inchi.startswith("InChI="):
                rdmol = Chem.MolFromInchi(smiles_or_inchi)
                if rdmol is not None:
                    yield str(mol_id), rdmol, Chem.MolToSmiles(rdmol)
            else:
                rdmol = Chem.MolFromSmiles(smiles_or_inchi)
                if rdmol is not None:
                    yield str(mol_id), rdmol, smiles_or_inchi


def sdf_mol_supplier(filename: str, gen_ids: bool=True, mol_id_prop: Union[str, None]=None) -> IterableType[Tuple[str, Chem.Mol, str]]:
    """Generator function that reads from a .sdf file.

    Parameters
    ----------
    filename : str
        .sdf filename.

    gen_ids: bool
        if True, automatically generate ids otherwise use from input file.

    mol_id_prop: str or None
        mol_id property in sdf file to use as id. if None, automatically generate ids using index in file.

    Yields
    -------
    tuple
        int id and rdkit mol.
    """
    if filename.endswith(".gz"):
        import gzip

        gzf = gzip.open(filename)
        suppl = Chem.ForwardSDMolSupplier(gzf)
    else:
        suppl = Chem.ForwardSDMolSupplier(filename)
    for new_mol_id, rdmol in enumerate(suppl, 1):
        if rdmol:
            mol_id = None
            if mol_id_prop is not None:
                try:
                    mol_id = rdmol.GetProp(mol_id_prop)
                except KeyError:
                    pass
                if mol_id is None:
                    print(f"WARNING: {mol_id_prop} not found in the sdf file! using the index as mol_id instead.")
                    
            if gen_ids or mol_id is None:
                mol_id = new_mol_id

            yield str(mol_id), rdmol, Chem.MolToSmiles(rdmol)
        else:
            continue


def get_mol_supplier(io_source: Any) -> Union[Callable[..., IterableType[Tuple[str, Chem.Mol]]], None]:
    """Returns a mol supplier depending on the object type and file extension.

    Parameters
    ----------
    mols_source : str or iterable
        .smi or .sdf filename or iterable.

    Returns
    -------
    callable
        function that will read the molecules from the input.
    """
    supplier = None
    if isinstance(io_source, str):
        split_source = io_source.split(".")
        if split_source[-1] == "gz":
            input_type = split_source[-2]
        else:
            input_type = split_source[-1]
        if input_type == "smi" or input_type == "inchi":
            supplier = smi_mol_supplier
        elif input_type == "sdf":
            supplier = sdf_mol_supplier
    elif isinstance(io_source, Iterable):
        supplier = it_mol_supplier
    else:
        raise Exception("Invalid input")
    return supplier
