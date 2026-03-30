"""
It contains a function to create a database from molecular fingerprints.

Author: Rajendra Kumar
"""

from typing import Any, Dict, List, Tuple
from datetime import datetime
import json
import sys
import glob
import os

from . import mol_fp_supplier, mol_fp_parallel_supplier, FPCalculator
from . import uffpsimLib

def create_database(input_file: str, db_file: str, fp_type: str, fp_params: Dict[str, Any] = None, 
                    gen_ids: bool = True, mol_id_prop=None, mol_id_max_chars: int = 15, info: dict[str, Any] | str = "",
                    inner_clustering_threshold: float = 0.2, cluster_mode: str = "memory", cluster_parallel: bool = False) -> None:
    """Creates a new database from the input file containing molecular fingerprints.

    Parameters
    ----------
    input_file : str
        The path to the input file containing molecular data.
    db_file : str
        The path to the database file where the fingerprints will be stored.
    fp_type : str
        The type of molecular fingerprint to be calculated.
    fp_params : Dict[str, Any], optional
        Additional parameters for fingerprint calculation.
    gen_ids : bool, optional
        Whether to generate unique molecule IDs if not provided in the input file.
    mol_id_prop : str, optional
        The property name to use as molecule ID from the input file.
    mol_id_max_chars : int, optional
        The maximum number of characters allowed for molecule IDs.
    info : dict[str, Any] | str, optional
        Additional information to be stored in the database.
    inner_clustering_threshold : float, optional
        The threshold for inner clustering.
    cluster_mode : str, optional
        The mode for storing and clustering the fingerprints. Currently, accepted keywords are `memory` and `disk`.
    cluster_parallel : bool, optional
        (Default: `False`). If `True`, enable clustering in parallel using OpenMP. The number of threads 
        can be controlled by setting `OMP_NUM_THRAEDS` variable.

    Returns
    -------
    None

    Raises
    ------
    Exception
        If any error occurs during the creation of the database.

    Example
    -------
    >>> # Create a database from a SDF file with Morgan fingerprints and inner clustering threshold 0.1
    >>> create_database(input_file='tests/data/10mols.sdf', db_file='from-sdf.h5', fp_type='Morgan', inner_clustering_threshold=0.1)
    >>> 
    >>> # Create a database from a SMILES file with Morgan fingerprints and inner clustering threshold 0.2
    >>> create_database(input_file='tests/data/10mols-with-ids.smi', db_file='from-smiles-with-ids.h5', fp_type='Morgan', gen_ids=False, inner_clustering_threshold=0.2)
    """

    fp_calculator = FPCalculator(fp_type, fp_params)
    fp_params_json = json.dumps({"fp_type": fp_calculator.type, "fp_params": fp_calculator.parameters})

    if isinstance(info, dict):
        info = json.dumps(info)

    fpstore = uffpsimLib.FingerprintStore(db_file, mol_id_max_chars=mol_id_max_chars, mode="w",
                                          fpSize=fp_calculator.parameters["fpSize"], fp_params=fp_params_json,
                                          info=info,cluster_threshold=inner_clustering_threshold,
                                          cluster_mode=cluster_mode,cluster_parallel=cluster_parallel)
    
    try:
        fps: List[Tuple[str, str]] = []
        total_processed = 0
        start_time = datetime.now()
        sys.stdout.write(" Started adding fingerprints...\n")
        for mol_id, fp, smiles in mol_fp_supplier(input_file, fp_type, fp_params=fp_params, 
                                                    gen_ids=gen_ids, mol_id_prop=mol_id_prop, 
                                                    mol_id_length=mol_id_max_chars):
            fps.append((mol_id, fp, smiles))

            # Appending to the fpstore object, it saves a lot of memory as string fp is converted to compact uint64 array 
            if len(fps)%100000 == 0:
                fpstore.append_fingerprints(fps)
                total_processed += len(fps)
                fps = []
                elapsed_time = datetime.now() - start_time
                sys.stdout.write(f"\r Added {total_processed} fingerprints. Total Time: {elapsed_time.total_seconds()/60:6.3f} mins...")
                sys.stdout.flush()

        if len(fps) > 0:
            total_processed += len(fps)
            fpstore.append_fingerprints(fps)
            elapsed_time = datetime.now() - start_time
        sys.stdout.write("\r Finished adding fingerprints\n")
        sys.stdout.write(f" Total fingerprints added: {total_processed};  Total Time: {elapsed_time.total_seconds()/60:6.3f} mins...\n")
        sys.stdout.flush()

        fpstore.perform_clustering_write()
    except Exception as e:
        sys.stderr.write(f"Error occurred while creating database: {str(e)}\n")
        fpstore.close()
        raise e

    fpstore.close()


def create_database_parallel(input_file: str, db_file: str, fp_type: str, workers= 4, fp_params: Dict[str, Any] = None, 
                    gen_ids: bool = True, mol_id_prop=None, mol_id_max_chars: int = 15, info: dict[str, Any] | str = "",
                    inner_clustering_threshold: float = 0.2, cluster_mode: str = "memory", cluster_parallel: bool = False) -> None:
    """Creates a new database from the input file containing molecular fingerprints using parallel processing.

    .. note:: It uses the `multiprocessing` module to create multiple worker processes. Therefore it does not scale well with increasing number of workers and requires more memory.

    Parameters
    ----------
    input_file : str
        The path to the input file containing molecular data.
    db_file : str
        The path to the database file where the fingerprints will be stored.
    fp_type : str
        The type of molecular fingerprint to be calculated.
    workers : int, optional
        The number of worker processes to use for parallel processing.
    fp_params : Dict[str, Any], optional
        Additional parameters for fingerprint calculation.
    gen_ids : bool, optional
        Whether to generate unique molecule IDs if not provided in the input file.
    mol_id_prop : str, optional
        The property name to use as molecule ID from the input file.
    mol_id_max_chars : int, optional
        The maximum number of characters allowed for molecule IDs.
    info : dict[str, Any] | str, optional
        Additional information to be stored in the database.
    inner_clustering_threshold : float, optional
        The threshold for inner clustering.
    cluster_mode : str, optional
        The mode for clustering the fingerprints. Currently, accepted keywords are `memory` and `disk`.
    cluster_parallel : bool
        (Default: `False`). If `True`, enable clustering in parallel using OpenMP. The number of threads 
        can be controlled by setting `OMP_NUM_THRAEDS` variable.    

    Returns
    -------
    None

    Raises
    ------
    Exception
        If any error occurs during the creation of the database.
    """
    fp_calculator = FPCalculator(fp_type, fp_params)
    fp_params_json = json.dumps({"fp_type": fp_calculator.type, "fp_params": fp_calculator.parameters})
    if isinstance(info, dict):
        info = json.dumps(info)

    fpstore = uffpsimLib.FingerprintStore(db_file, mol_id_max_chars=mol_id_max_chars, mode="w",
                                          fpSize=fp_calculator.parameters["fpSize"], fp_params=fp_params_json,
                                          cluster_threshold=inner_clustering_threshold,
                                          cluster_mode=cluster_mode,cluster_parallel=cluster_parallel)
    
    try:
        fps: List[Tuple[str, str]] = []
        total_processed = 0
        start_time = datetime.now()
        sys.stdout.write(" Started adding fingerprints...\n")
        for mol_id, fp, smiles in mol_fp_parallel_supplier(input_file, fp_type, workers=workers, fp_params=fp_params, 
                                                    gen_ids=gen_ids, mol_id_prop=mol_id_prop, 
                                                    mol_id_length=mol_id_max_chars):
            fps.append((mol_id, fp, smiles))

            # Appending to the fpstore object, it saves a lot of memory as string fp is converted to compact uint64 array 
            if len(fps)%100000 == 0:
                fpstore.append_fingerprints(fps)
                total_processed += len(fps)
                fps = []
                elapsed_time = datetime.now() - start_time
                sys.stdout.write(f"\r Added {total_processed} fingerprints. Total Time: {elapsed_time.total_seconds()/60:6.3f} mins...")
                sys.stdout.flush()

        if len(fps) > 0:
            total_processed += len(fps)
            fpstore.append_fingerprints(fps)
            elapsed_time = datetime.now() - start_time
        sys.stdout.write("\n Finished adding fingerprints\n")
        sys.stdout.write(f" Total fingerprints added: {total_processed};  Total Time: {elapsed_time.total_seconds()/60:6.3f} mins...\n")
        sys.stdout.flush()

        fpstore.perform_clustering_write()
    except Exception as e:
        sys.stderr.write(f"Error occurred while creating database: {str(e)}\n")
        fpstore.close()
        raise e

    fpstore.close()

def create_database_from_files(input_files: List[str], db_file: str, fp_type: str, workers= 1, fp_params: Dict[str, Any] = None, 
                    gen_ids: bool = True, mol_id_prop=None, mol_id_max_chars: int = 15, info: dict[str, Any] | str = "",
                    inner_clustering_threshold: float = 0.2, cluster_mode: str = "memory", cluster_parallel: bool = False):
    
    fp_calculator = FPCalculator(fp_type, fp_params)
    fp_params_json = json.dumps({"fp_type": fp_calculator.type, "fp_params": fp_calculator.parameters})
    if isinstance(info, dict):
        info = json.dumps(info)

    fpstore = uffpsimLib.FingerprintStore(db_file, mol_id_max_chars=mol_id_max_chars, mode="w",
                                          fpSize=fp_calculator.parameters["fpSize"], fp_params=fp_params_json,
                                          cluster_threshold=inner_clustering_threshold,
                                          cluster_mode=cluster_mode,cluster_parallel=cluster_parallel)
    
    sys.stdout.write(" Started adding fingerprints...\n")
    fps: List[Tuple[str, str]] = []
    total_processed = 0
    start_time = datetime.now()
    total_files = len(input_files)
    for n, input_file in enumerate(input_files, start=1):
        try:
            for mol_id, fp, smiles in mol_fp_parallel_supplier(input_file, fp_type, workers=workers, fp_params=fp_params, 
                                                        gen_ids=gen_ids, mol_id_prop=mol_id_prop, 
                                                        mol_id_length=mol_id_max_chars):
                fps.append((mol_id, fp, smiles))

                # Appending to the fpstore object, it saves a lot of memory as string fp is converted to compact uint64 array 
                if len(fps)%100000 == 0:
                    fpstore.append_fingerprints(fps)
                    total_processed += len(fps)
                    fps = []
                    elapsed_time = datetime.now() - start_time
                    sys.stdout.write(f"\r Processing file [{n}/{total_files}]: {input_file}, Added {total_processed} fingerprints. Total Time: {elapsed_time.total_seconds()/60:6.3f} mins...")
                    sys.stdout.flush()
            
        except Exception as e:
            sys.stderr.write(f"Error occurred while creating database: {str(e)}\n")
            fpstore.close()
            raise e

    try:
        # last of the fingerprints
        if len(fps) > 0:
            total_processed += len(fps)
            fpstore.append_fingerprints(fps)
            elapsed_time = datetime.now() - start_time
        sys.stdout.write("\n Finished adding fingerprints\n")
        sys.stdout.write(f" Total fingerprints added: {total_processed};  Total Time: {elapsed_time.total_seconds()/60:6.3f} mins...\n")
        sys.stdout.flush()

        # now clustering
        fpstore.perform_clustering_write()

    except Exception as e:
        sys.stderr.write(f"Error occurred while performing inner clustering: {str(e)}\n")
        fpstore.close()
        raise e
    
    fpstore.close()

def create_database_from_dir(input_dir: str, db_file: str, fp_type: str, suffix = "sdf.gz", workers= 1, fp_params: Dict[str, Any] = None, 
                    gen_ids: bool = True, mol_id_prop=None, mol_id_max_chars: int = 15, info: dict[str, Any] | str = "",
                    inner_clustering_threshold: float = 0.2, cluster_mode: str = "memory", cluster_parallel: bool = False):
    

    if suffix not in ["sdf", "sdf.gz", "smi"]:
        raise ValueError(f"Unsupported file suffix: {suffix}. Supported suffixes are: 'sdf', 'sdf.gz', 'smi'.")
    
    input_files = glob.glob(os.path.join(input_dir, f'*.{suffix}'))
    create_database_from_files(input_files, db_file, fp_type, workers=workers, fp_params=fp_params, 
                    gen_ids=gen_ids, mol_id_prop=mol_id_prop, mol_id_max_chars=mol_id_max_chars,
                    info=info, inner_clustering_threshold=inner_clustering_threshold,
                    cluster_mode=cluster_mode, cluster_parallel=cluster_parallel)

def redo_inner_clustering(db_file: str, threshold: float, cluster_mode: str = "memory", cluster_parallel: bool = False) -> None:
    """Redoes the inner clustering of the database using the specified threshold.

    Parameters
    ----------
    db_file : str
        The path to the database file.
    threshold : float
        The threshold for inner clustering.
    cluster_mode : str, optional
        The mode for clustering the fingerprints. Currently, accepted keywords are `memory` and `disk`.
    cluster_parallel : bool
        (Default: `False`). If `True`, enable clustering in parallel using OpenMP. The number of threads 
        can be controlled by setting `OMP_NUM_THRAEDS` variable.

    Returns
    -------
    None

    Raises
    ------
    Exception
        If any error occurs during the redoing of inner clustering.
    """
    fpstore = uffpsimLib.FingerprintStore(db_file, mode="a", cluster_mode = cluster_mode, cluster_parallel= cluster_parallel)
    fpstore.redo_clustering_write(threshold)
    del fpstore


def build_mol_id_index_table(db_file: str) -> None:
    """Builds the serialized MolIdIndexTable for an existing uffpsim HDF5 database.

    Parameters
    ----------
    db_file : str
        The path to the database file.

    Returns
    -------
    None

    Raises
    ------
    Exception
        If any error occurs while building the serialized molecule ID index table.
    """
    fpstore = uffpsimLib.FingerprintStore(db_file, mode="a")
    try:
        fpstore.build_mol_id_index_table()
    except Exception as e:
        fpstore.close()
        raise e

    fpstore.close()

def update_database(input_file: str, db_file: str, mol_id_prop=None, info: dict[str, Any] | None=None):
    """Updates the database with new molecular fingerprints from the input file.

    It checks for new molecule IDs in the input file and updates the database accordingly. If
    molecule IDs are present in the database, the existing fingerprints are updated with the new one.

    .. note:: it also performs inner clustering on the updated pop-count bins.

    Parameters
    ----------
    input_file : str
        The path to the input file containing molecular data.
    db_file : str
        The path to the database file.
    mol_id_prop : str, optional
        The property name to use as molecule ID from the input file.
    info : dict[str, Any] | None, optional
        Additional information to be stored in the database.

    Returns
    -------
    None

    Raises
    ------
    Exception
        If any error occurs while updating the database.
    """
    fpstore = uffpsimLib.FingerprintStore(db_file, mode="a")
    fp_input_arguments = json.loads(fpstore.fp_params_json)
    fp_type = fp_input_arguments["fp_type"]
    fp_params = fp_input_arguments["fp_params"]

    try:
        input_fps_data = []
        for mol_id, fp, smiles in mol_fp_supplier(input_file, fp_type, fp_params=fp_params, 
                                                    gen_ids=False, mol_id_prop=mol_id_prop, 
                                                    mol_id_length=fpstore.mol_id_max_chars):
            input_fps_data.append((mol_id, fp, smiles))

        fpstore.update_fingerprints(input_fps_data)

        if info is not None:
            old_info = json.loads(fpstore.info)
            new_info = {**old_info, **info}
            fpstore.set_info(json.dumps(new_info))

    except Exception as e:
        fpstore.close()
        raise e
    
    fpstore.close()

def get_database_info(db_file: str) -> Dict[str, Any]:
    """Retrieves the information stored in the database.

    Parameters
    ----------
    db_file : str
        The path to the database file.

    Returns
    -------
    Dict[str, Any]
        The information stored in the database.

    Raises
    ------
    Exception
        If any error occurs while retrieving the information from the database.
    """
    fpstore = uffpsimLib.FingerprintStore(db_file, mode="r")
    info = {}
    try:
        info = json.loads(fpstore.info)
    except:
        pass
    fpstore.close()
    return info