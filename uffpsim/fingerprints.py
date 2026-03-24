"""
Calculate molecular fingerprints.
Most part taken from FPSim2 source code.
"""

from typing import Any, Dict
from collections.abc import Iterable
import re
from rdkit import rdBase
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import rdFingerprintGenerator
from rdkit.Avalon import pyAvalonTools
from rdkit import Chem

blocker = rdBase.BlockLogs()

MOLFILE_RE = r" [vV][23]000$"

FP_TYPES = [ "Morgan", "RDKit", "AtomPair", "TopologicalTorsion", "MACCSKeys", "Avalon", "RDKPatternFingerprint"]

FP_TYPE_GENERATORS = {
    "Morgan": rdFingerprintGenerator.GetMorganGenerator,
    "RDKit": rdFingerprintGenerator.GetRDKitFPGenerator,
    "AtomPair": rdFingerprintGenerator.GetAtomPairGenerator,
    "TopologicalTorsion": rdFingerprintGenerator.GetTopologicalTorsionGenerator,
}

FP_TYPE_FUNCS = {
    "MACCSKeys": rdMolDescriptors.GetMACCSKeysFingerprint,
    "Avalon": pyAvalonTools.GetAvalonFP,
    "RDKPatternFingerprint": Chem.PatternFingerprint,
}

FP_TYPE_DEFAULT_PARAMETERS = {
    "MACCSKeys": {},
    "Avalon": {
        "nBits": 512,
        "isQuery": False,
        "resetVect": False,
        "bitFlags": 15761407,
    },
    "Morgan": {
        "radius": 2,
        "fpSize": 2048, 
        "countSimulation": False,
        "includeChirality": False,
        "useBondTypes": True,
        "onlyNonzeroInvariants": False,
        "includeRingMembership": True,
        "countBounds": None,
        "atomInvariantsGenerator": None, 
        "bondInvariantsGenerator": None, 
        "includeRedundantEnvironments": False
    },
    "TopologicalTorsion": {
        "fpSize": 2048,
        "torsionAtomCount": 4,
        "countSimulation": True,
        "countBounds": None,
        "atomInvariantsGenerator": None,
        "includeChirality": False,
    },
    "AtomPair": {
        "fpSize": 2048,
        "minDistance": 1,
        "maxDistance": 30,
        "includeChirality": False,
        "use2D": True,
        "countSimulation": True,
        "countBounds": None,
        "atomInvariantsGenerator": None
    },
    "RDKit": {
        "minPath": 1,
        "maxPath": 7,
        "fpSize": 2048,
        "branchedPaths": True,
        "useBondOrder": True,
        "countSimulation": False,
        "countBounds": None,
        "numBitsPerFeature": 2,
        "atomInvariantsGenerator": None
    },
    "RDKPatternFingerprint": {"fpSize": 2048, "atomCounts": [], "setOnlyBits": None},
}

class FPCalculator:
    """A class for calculating fingerprints."""
    def __init__(self, fp_type: str, fp_params: Dict[str, Any] = None) -> None:
        self.type = None
        self.parameters = None

        self._set_fp_type(fp_type)
        self._set_fp_parameters(fp_params)

        self.fp_generator = None
        self.fp_function = None
        self.calculator = None
        if self.type in FP_TYPE_GENERATORS:
            self.fp_generator = FP_TYPE_GENERATORS[self.type](**self.parameters)
            self.calculator = lambda mol: self.fp_generator.GetFingerprint(mol).ToBitString()
        else:
            self.fp_function = FP_TYPE_FUNCS[self.type]
            self.calculator = lambda mol: self.fp_function(mol, **self.parameters).ToBitString()

        self.fpSize = None
        if "fpSize" in self.parameters:
            self.fpSize = self.parameters["fpSize"]
        if "nBits" in self.parameters:
            self.fpSize = self.parameters["nBits"]

        if self.fpSize is not None:
            if self.fpSize % 64 != 0:
                raise ValueError(f"Number of bits in the fingerprint must be multiple of 64 to exploit the full performance. Here it is {self.fpSize}...")
        else:
            raise ValueError("Cannot determine number of bits in the Fingerprint.")


    def _set_fp_type(self, fp_type) -> None:
        if not (fp_type in FP_TYPE_FUNCS or fp_type in FP_TYPE_GENERATORS):
            raise ValueError("Invalid fingerprint type! Please choose one of the following: " + ", ".join(FP_TYPE_FUNCS.keys()))
        self.type = fp_type

    def _set_fp_parameters(self, fp_params: Dict[str, Any] = None) -> None:
        if not fp_params:
            self.parameters = FP_TYPE_DEFAULT_PARAMETERS[self.type]
        else:
            for key in fp_params:
                if key not in FP_TYPE_DEFAULT_PARAMETERS[self.type]:
                    raise ValueError("Invalid parameter for fingerprint type! Please choose one of the following: " + 
                                     ", ".join(FP_TYPE_DEFAULT_PARAMETERS[self.type].keys()))
                
            self.parameters = {**FP_TYPE_DEFAULT_PARAMETERS[self.type], **fp_params}

    def __call__(self, mol: Chem.rdchem.Mol) -> str:
        return self.calculator(mol)

def get_fp_length(fp_type: str, fp_params: Dict[str, Any]) -> int:
    """Returns the FP length given the name of the FP function and it's parameters.

    Parameters
    ----------
    fp_type : str
         Name of the function used to generate the fingerprints.

    fp_params: dict
        Parameters used to generate the fingerprints.

    Returns
    -------
    fp_length: int
        fp length of the fingerprint.
    """
    fp_length = None
    if "nBits" in fp_params.keys():
        fp_length = fp_params["nBits"]
    elif "fpSize" in fp_params.keys():
        fp_length = fp_params["fpSize"]
    if fp_type == "MACCSKeys":
        fp_length = 166
    if not fp_length:
        raise Exception("fingerprint size is not specified")
    return fp_length

def load_molecule(molecule: Any) -> tuple[Chem.Mol, str]:
    """Reads SMILES, molblock or InChI and returns a RDKit mol.

    Parameters
    ----------
    molecule : Any
         Chem.Mol, SMILES, molblock or InChI.

    Returns
    -------
    mol: ROMol
        RDKit molecule.
    """
    smiles = ''
    if isinstance(molecule, Chem.Mol):
        smiles = Chem.MolToSmiles(molecule)
        return (molecule, smiles)
    if re.search(MOLFILE_RE, molecule, flags=re.MULTILINE):
        rdmol = Chem.MolFromMolBlock(molecule)
        return (rdmol, Chem.MolToSmiles(rdmol))
    elif molecule.startswith("InChI="):
        rdmol = Chem.MolFromInchi(molecule)
        return (rdmol, Chem.MolToSmiles(rdmol))
    else:
        rdmol = Chem.MolFromSmiles(molecule)
        return (rdmol, f'{molecule}')