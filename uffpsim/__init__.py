import os
from . import uffpsimLib
from .fingerprints import FPCalculator, FP_TYPE_DEFAULT_PARAMETERS
from .fp_supplier import mol_fp_supplier, it_mol_supplier, smi_mol_supplier, sdf_mol_supplier, get_mol_supplier
from .fp_parallel_supplier import mol_fp_parallel_supplier
from .database import create_database, create_database_parallel, update_database, redo_inner_clustering, get_database_info
from .search_engine import UFFPSimSearchEngine
import importlib.metadata

os.environ["OMP_CANCELLATION"] = "true"

__all__ = ["uffpsimLib", "mol_fp_supplier", "FPCalculator", "FP_TYPE_DEFAULT_PARAMETERS", "mol_fp_parallel_supplier",
           "create_database", "create_database_parallel", "update_database", "get_database_info",
            "redo_inner_clustering" ,"UFFPSimSearchEngine"]
__version__ = importlib.metadata.version(__package__ or __name__)
