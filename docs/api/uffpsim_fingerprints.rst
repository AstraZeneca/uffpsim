Fingerprints Module
===================

.. currentmodule:: uffpsim.fingerprints

Module Overview
---------------

The fingerprints module provides fingerprint calculation support for multiple fingerprint types, including Morgan, RDKit, AtomPair, TopologicalTorsion, MACCSKeys, Avalon, and RDKPatternFingerprint. It also includes utility functions for loading molecules from SMILES, InChI, or molblock formats.

Classes
-------

.. autoclass:: FPCalculator
   :members: __init__, __call__
   :undoc-members:
   :show-inheritance:

Constants
---------

.. autodata:: FP_TYPE_DEFAULT_PARAMETERS
   :no-value:

.. autodata:: FP_TYPES
   :no-value:

Functions
---------

.. autofunction:: get_fp_length

.. autofunction:: load_molecule
