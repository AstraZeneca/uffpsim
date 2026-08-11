Demo Web-App for pre-built Chembl and PubChem Databases
=======================================================

This page provides guides for running uffpsim's built-in web application for pre-built Chembl and PubChem databases.

Prerequisites
-------------

Before starting, install uffpsim and its dependencies:

.. code-block:: bash

    UFFPSIM_NATIVE=1 pip install -v ".[web]"

Chembl Database (~2.6 Million Compounds)
-------------------------------------------

Step 1: Download a Database File
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The Chembl database file can be downloaded from this `link <https://huggingface.co/datasets/rjdkmr/UFFPSim/resolve/main/chembl/database/chembl35_b2048_ct16.h5?download=true>`_ or as follows:

.. code-block:: bash

    curl -L https://huggingface.co/datasets/rjdkmr/UFFPSim/resolve/main/chembl/database/chembl35_b2048_ct16.h5?download=true -o ./chembl_35_b2048.h5

Step 2: Launch the Web App
~~~~~~~~~~~~~~~~~~~~~~~~~~

Run the following command from the uffpsim repository root:

.. code-block:: bash

    uffpsim launch-web-app -d chembl_35_b2048.h5


Step 3: Access the Web Application
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Open a browser and navigate to:

    http://127.0.0.1:5000/

You can now enter SMILES strings in the search form and view similarity results with molecular structure images.


PubChem Database (~124 Million Compounds)
--------------------------------------------------

Step 1: Download a Database File
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The PubChem database file can be downloaded from this `link <https://huggingface.co/datasets/rjdkmr/UFFPSim/resolve/main/pubchem/pubchem_b512_0_27.h5?download=true>`_ or as follows:

.. code-block:: bash

    curl -L https://huggingface.co/datasets/rjdkmr/UFFPSim/resolve/main/pubchem/pubchem_b512_0_27.h5?download=true -o ./pubchem_b512.h5


Step 2: Launch the Web App with Full Results Mode
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Use ``--results full`` to return compound IDs, SMILES strings, and molecular images for each hit.
Build the molecule ID index first so SMILES lookup is available:


.. code-block:: bash

    uffpsim launch-web-app -d pubchem_b512.h5 --results full

Step 3: Access the Web Application
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Navigate to:

    http://127.0.0.1:5000/

The web app will display full hit details including molecular structures for each matching compound.
