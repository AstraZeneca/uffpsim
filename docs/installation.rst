Installation
============

The simplest way to install uffpsim is:

.. code-block:: bash

    UFFPSIM_NATIVE=1 pip install -v .

If the installation fails due to missing dependencies or shared libraries, you can alternatively use a conda environment. The following will create a new conda environment in the current working directory and install all requirements from the ``dev-environment.yaml`` file.

.. code-block:: bash

    # Create development conda environment
    conda env create --prefix ./venv --file dev-environment.yaml

    # Activate the environment
    conda activate ./venv

    # Install uffpsim
    UFFPSIM_NATIVE=1 pip install -v .

Installation with ``-march=native`` support
--------------------------------------------

This will compile and build ``uffpsim`` that will be highly optimized for the current CPU.

.. code-block:: bash

    UFFPSIM_NATIVE=1 pip install -v .

Installation with AVX512-VPOPCNTDQ support
--------------------------------------------

The `AVX512-VPOPCNTDQ <https://en.wikipedia.org/wiki/AVX-512#VPOPCNTDQ_and_BITALG>`_ instruction set could speed-up similarity calculation by more than 50% on newest Intel CPUs.

To check whether CPU supports AVX512-VPOPCNTDQ, the following command on Linux should show ``avx512_vpopcntdq`` in the list.

.. code-block:: bash

    lscpu | grep Flags | grep avx

If ``avx512_vpopcntdq`` is not in the above list, it means the CPU does not have the capability to use this method.

This will compile and build ``uffpsim`` with AVX512-VPOPCNTDQ support.

.. code-block:: bash

    UFFPSIM_AVX512=1 pip install -v .

Development Setup
-----------------

For development, run:

.. code-block:: bash

    conda env create --prefix ./venv --file dev-environment.yaml  # Create development conda environment
    conda activate ./venv                                           # Activate the environment
    UFFPSIM_NATIVE=1 pip install -ve ".[web]"                       # Install uffpsim in editable mode
