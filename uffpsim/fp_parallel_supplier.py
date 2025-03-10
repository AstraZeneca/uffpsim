"""
Fingerprint suppliers in parallel processing

Author: Rajendra Kumar
"""

from typing import Any, Callable, Iterable as IterableType, Dict, List, Tuple, Union
from collections.abc import Iterable
from rdkit import rdBase
from rdkit import Chem
import io
from concurrent.futures import ProcessPoolExecutor
import uuid

from .fingerprints import FPCalculator, load_molecule

GlobalFPCalculators = {}

def sdf_supplier(io_resource):
    """
    A supplier function that reads a SDF (Structural Data File) file and yields blocks of data.

    Parameters
    ----------
    io_resource : str
        The path to the SDF file. The file can be compressed with gzip.

    Yields
    ------
    str
        A block of data from the SDF file.

    Raises
    ------
    Exception
        If the file cannot be opened or read.
    """
    if io_resource.endswith(".gz"):
        import gzip
        reader = gzip.open(io_resource, 'rt')
    else:
        reader = open(io_resource)

    block = []
    for line in reader:
        if "$$$$" in line:
            block.append(line)
            yield ''.join(block)

            block = []
        else:
            block.append(line)

    if (len(block) > 0):
        yield ''.join(block)

def smi_supplier(io_resource):
    """A supplier function that reads a SMILES (Simplified Molecular Input Line Entry) file and yields blocks of data.

    Parameters
    ----------
    io_resource : str
        The path to the SMILES file. The file can be compressed with gzip.

    Yields
    ------
    str
        A block of data from the SMILES file. Each block contains one SMILES string.

    Raises
    ------
    Exception
        If the file cannot be opened or read.
    """
    with open(io_resource, "r") as f:
        for line in f:
            yield line.strip()

def iterable_supplier(io_resource):
    """A supplier function that reads an iterable object and yields blocks of data.

    Parameters
    ----------
    io_resource : iterable
        An iterable object. Each element in the iterable is considered as a block of data.

    Yields
    ------
    Any
        A block of data from the iterable.

    Raises
    ------
    Exception
        If the iterable object cannot be iterated over.
    """
    for block in io_resource:
        yield block

def supplier_factory(io_resource):
    """
    A factory function that creates a supplier based on the type of the input resource.

    Parameters
    ----------
    io_resource : Any
        The input resource. It can be a file path (str), an iterable object, or any other type.

    Returns
    -------
    Tuple[str, Callable]
        A tuple containing the type of the input resource (io_type) and the corresponding supplier function.

    Raises
    ------
    Exception
        If the input resource type is not supported.
    """
    supplier = None
    io_type = None
    if isinstance(io_resource, str):
        split_source = io_resource.split(".")
        if split_source[-1] == "gz":
            input_type = split_source[-2]
        else:
            input_type = split_source[-1]

        if input_type == "smi" or input_type == "inchi":
            supplier = smi_supplier
            io_type = "smi"
        elif input_type == "sdf":
            supplier = sdf_supplier
            io_type = "sdf"
    elif isinstance(io_resource, Iterable):
        supplier = iterable_supplier
        io_type = "iterable"
    else:
        raise Exception("Invalid input")
    return io_type, supplier

def process_smi_block(index, block, config):
    """A function that processes a block of SMILES strings and generates molecular identifiers, RDKit molecules, and SMILES strings.

    Parameters
    ----------
    index : int
        The index of the block.
    block : str
        The block of SMILES strings.
    config : dict
        A dictionary containing configuration parameters.

    Returns
    -------
    Tuple[str, Chem.Mol, str]
        A tuple containing the molecular identifier (mol_id), the corresponding RDKit molecule (rdmol), and the SMILES string.

    Raises
    ------
    ValueError
        If the mol_id is too long.
    """
    gen_ids = config["gen_ids"]
    mol = block.split()
    if  gen_ids:
        smiles_or_inchi = mol[0]
        mol_id = index
    else:
        if len(mol) == 1:
            smiles_or_inchi = mol[0]
            mol_id = index
            if config["show_warning"]:
                config["show_warning"] = False
                print("WARNING: No mol_id found in the SMILES file!"
                                "Consider setting gen_ids=True when running to automatically generate the mol_id from index.")
        else:
            smiles_or_inchi = mol[0]
            mol_id = mol[1]
    
    smiles_or_inchi = smiles_or_inchi.strip()
    if smiles_or_inchi.startswith("InChI="):
        rdmol = Chem.MolFromInchi(smiles_or_inchi)
        if rdmol is None:
            return str(mol_id), None, None
        else:
            return str(mol_id), rdmol, Chem.MolToSmiles(rdmol)
    else:
        rdmol = Chem.MolFromSmiles(smiles_or_inchi)
        return str(mol_id), rdmol, smiles_or_inchi

def process_sdf_block(index, block, config):
    """
    A function that processes a block of SDF strings and generates molecular identifiers, RDKit molecules, and SMILES strings.

    Parameters
    ----------
    index : int
        The index of the block.
    block : str
        The block of SDF strings.
    config : dict
        A dictionary containing configuration parameters.

    Returns
    -------
    Tuple[str, Chem.Mol, str]
        A tuple containing the molecular identifier (mol_id), the corresponding RDKit molecule (rdmol), and the SMILES string.

    Raises
    ------
    ValueError
        If the mol_id is too long.
    """
    mol_id_prop = config["mol_id_prop"]
    gen_ids = config["gen_ids"]
    suppl = Chem.ForwardSDMolSupplier(io.BytesIO(block.encode('utf-8')))
    mol = list(suppl)[0]

    if mol is None:
        return str(index), None, None
    
    mol_id = None
    if mol_id_prop is not None:
        try:
            mol_id = mol.GetProp(mol_id_prop)
        except:
            pass
        if mol_id is None:
            if config["show_warning"]:
                config["show_warning"] = False
                print(f"WARNING: {mol_id_prop} not found in the sdf file! using the index as mol_id instead.")
            
    if gen_ids or mol_id is None:
        mol_id = index

    return str(mol_id), mol, Chem.MolToSmiles(mol)

def process_iterable_entity(index, entity, config):
    """A function that processes a block of iterable entities and generates molecular identifiers, RDKit molecules, and SMILES strings.

    Parameters
    ----------
    index : int
        The index of the block.
    entity : iterable
        The block of iterable entities. Each entity can be a string or a tuple containing a molecule string and its identifier.
    config : dict
        A dictionary containing configuration parameters.

    Returns
    -------
    Tuple[str, Chem.Mol, str]
        A tuple containing the molecular identifier (mol_id), the corresponding RDKit molecule (rdmol), and the SMILES string.

    Raises
    ------
    ValueError
        If the mol_id is too long.
    """
    if isinstance(entity, str):
        mol_string = entity
        mol_id = index
    else:
        if config["gen_ids"]:
            mol_string = entity[0]
            mol_id = index
        else:
            if len(entity) == 1:
                mol_string = entity
                mol_id = index
                if config["show_warning"]:
                    config["show_warning"] = False
                    raise  Warning("No mol_id found in the iterable!"
                                    "Consider setting gen_ids=True when running to automatically generate the mol_id from index.")
                    
            else:
                mol_string = entity[0]
                mol_id = entity[1]
        
    (rdmol, smiles) = load_molecule(mol_string)
    return str(mol_id), rdmol, smiles

def process_block(data):
    """
    A function that processes a block of data and generates molecular identifiers, RDKit molecules, and SMILES strings.

    Parameters
    ----------
    data : tuple
        A tuple containing the index blocks, the type of input data (io_type), the UUID of the current process, and the configuration parameters.

    Returns
    -------
    List[Tuple[str, Chem.Mol, str]]
        A list of tuples, where each tuple contains the molecular identifier (mol_id), the corresponding RDKit molecule (rdmol), and the SMILES string.

    Raises
    ------
    ValueError
        If the mol_id is too long.
    """
    blocker = rdBase.BlockLogs()
    (index_blocks, io_type, uuid_value, config) = data
    fp_calculator = GlobalFPCalculators[uuid_value]
    process_func = None
    if io_type == "smi":
        process_func =  process_smi_block
    elif io_type == "sdf":
        process_func = process_sdf_block
    elif io_type == "iterable":
        process_func = process_iterable_entity
    else:
        raise ValueError("Unsupported format")
    
    out = []
    for index, block in index_blocks:
        mol_id, rdmol, smiles = process_func(index, block, config)
        if len(mol_id) > config["mol_id_length"]:
            raise ValueError(f"mol_id {mol_id} is too long. Please choose a smaller mol_id_length.")
    
        if rdmol is not None:
            out.append((mol_id, fp_calculator(rdmol), smiles))

    return out
        
def get_chunk_size(length, workers):
    """
    A function that calculates the optimal chunk size for processing data in parallel.

    Parameters
    ----------
    length : int
        The total number of data items to be processed.
    workers : int
        The number of worker processes available for parallel processing.

    Returns
    -------
    int
        The optimal chunk size for parallel processing.

    """
    if (length % workers == 0):
        chunks = int(length/workers)
    else:
        chunks = int(length/workers) + workers
    return chunks

def run_process(index_blocks, chunks, executor, io_type, curr_uuid, config):
    """
    A function that processes a list of index blocks in parallel using a ProcessPoolExecutor.

    Parameters
    ----------
    index_blocks : list
        A list of tuples, where each tuple contains the index and the corresponding block of data.
    chunks : int
        The size of each chunk for parallel processing.
    executor : concurrent.futures.ProcessPoolExecutor
        An instance of ProcessPoolExecutor for parallel processing.
    io_type : str
        The type of input data (e.g., "smi", "sdf", "iterable").
    curr_uuid : uuid.UUID
        The UUID of the current process.
    config : dict
        A dictionary containing configuration parameters.

    Returns
    -------
    list
        A list of tuples, where each tuple contains the molecular identifier (mol_id), the corresponding fingerprint (fp), and the SMILES string.

    """
    data= []
    for i in range(0, len(index_blocks), chunks):
        start = i
        end = i + chunks
        data.append((index_blocks[start:end], io_type, curr_uuid, config))

    out = []
    for out_i in executor.map(process_block, data):
        out.extend(out_i)
    return out

def mol_fp_parallel_supplier(input_file: Any, fp_type: str, workers: int = 1, gen_ids: bool=True, mol_id_prop: Union[str, None]=None, fp_params: Dict[str, Any] = None, mol_id_length: int = 15):
    """
    A function that reads a file containing molecular data (e.g., SMILES, SDF), generates molecular identifiers, calculates fingerprints, and yields the results in parallel using multiple worker processes.

    Parameters
    ----------
    input_file : Any
        The input file containing molecular data. It can be a file path (str), an iterable object, or any other type.
    fp_type : str
        The type of fingerprint to be calculated (e.g., "ECFP4", "MorganFP").
    workers : int, optional
        The number of worker processes to be used for parallel processing. Default is 1.
    gen_ids : bool, optional
        A flag indicating whether to generate molecular identifiers if they are not provided in the input file. Default is True.
    mol_id_prop : Union[str, None], optional
        The property name in the input file to be used as the molecular identifier. If None, the index of the molecule in the file will be used as the identifier. Default is None.
    fp_params : Dict[str, Any], optional
        Additional parameters for fingerprint calculation. Default is an empty dictionary.
    mol_id_length : int, optional
        The maximum length of the molecular identifier. If the generated identifier exceeds this length, a ValueError will be raised. Default is 15.

    Yields
    ------
    Tuple[str, Any, str]
        A tuple containing the molecular identifier (mol_id), the corresponding fingerprint (fp), and the SMILES string.

    Raises
    ------
    ValueError
        If the mol_id_length is not a multiple of 8 or if the input file cannot be opened or read.
    """
    config = {"gen_ids": gen_ids, "mol_id_prop": mol_id_prop, "show_warning": True, "mol_id_length": mol_id_length}
    
    if (mol_id_length+1) % 8 != 0:
        raise ValueError(f"mol_id_length {mol_id_length} must be less than one multiple of 8. e.g. 7, 15, 23, 31 etc.")
    
    fp_calculator = FPCalculator(fp_type, fp_params)
    
    curr_uuid = uuid.uuid4()
    chunksize = get_chunk_size(100000, workers)

    with ProcessPoolExecutor(max_workers=workers) as executor:
        try:
            GlobalFPCalculators[curr_uuid] = fp_calculator
            index_blocks = []
            io_type, supplier = supplier_factory(input_file)
            for index, block in enumerate(supplier(input_file)):
                index_blocks.append((index + 1, block))

                if len(index_blocks) >= 100000:
                    for mol_id, fp, smiles in run_process(index_blocks, chunksize, executor, io_type, curr_uuid, config):
                        yield mol_id, fp, smiles
                    index_blocks = []

            for mol_id, fp, smiles in run_process(index_blocks, chunksize, executor, io_type, curr_uuid, config):
                yield mol_id, fp, smiles

        except Exception as e:
            del GlobalFPCalculators[curr_uuid]
            raise e