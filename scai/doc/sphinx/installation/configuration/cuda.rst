Using CUDA in LAMA
^^^^^^^^^^^^^^^^^^

CUDA is needed to utilize CUDA capable GPUs from NVidia. The CUDA toolkit can be downloaded `here`__.

__ https://developer.nvidia.com/cuda-downloads

CUDA is optional and LAMA can be built without it. But you will not be able to take advantage of GPUs.

The configuration of LAMA usually finds an available CUDA installation on your system.
If not, you can give it a hint where to find it

.. code-block::

   cmake -D CUDA_TOOLKIT_ROOT=<path/to/cuda/installation>

If CUDA is available on your system but you do not want to use it, you can switch off its use as follows

.. code-block:: bash

   cmake -D USE_CUDA=OFF

Furthermore, you can change relevant CMake variables for CUDA by using the ccmake utility.

Beside the CUDA compiler, LAMA uses also:

- Thrust, a C++ template library for CUDA based on the Standard Template Library (STL). 
  LAMA uses it to implement some operations on matrices and vectors with minimal programming effort
  through a high-level interface that is fully interoperable with CUDA C.
  Since Thrust is a template library of header files, no further installation is necessary for using Thrust.

- the CUBLAS library

- the CUSparse library

Known problems:

Thrust might be confused about system files

.. code-block:: bash

   error: kernel launches from templates are not allowed in system files

Workaround

.. code-block:: bash

   unset CPLUS_INCLUDE_PATH
