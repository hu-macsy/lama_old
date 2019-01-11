.. _ParMetisPartitioning:

ParMetis Partitioning
^^^^^^^^^^^^^^^^^^^^^

|PARMETIS| is an MPI-based parallel library that implements a variety of algorithms for partitioning unstructured graphs, meshes, 
and for computing fill-reducing orderings of sparse matrices. ParMETIS extends the functionality provided by METIS and 
includes routines that are especially suited for parallel AMR computations and large scale numerical simulations. 

.. |PARMETIS| raw:: html

  <a href="http://glaros.dtc.umn.edu/gkhome/metis/parmetis/overview" target="_blank">ParMetis</a>

If your ParMetis installation is not in the default path define the PARMETIS_ROOT in the cmake call with -DPARMETIS_ROOT=<path_to_parmetis> 
or define the environment variables PARMETIS_INCLUDE_PATH and PARMETIS_LIBRARY_PATH as follows:

.. code-block:: bash

   export PARMETIS_INCLUDE_PATH=<path_to_parmetis>/include
   export PARMETIS_LIBRARY_PATH=<path_to_parmetis>/lib

Download latest stable PARMETIS release (e.g. parmetis-4.0.3.tar.gz), unpack it.

Before compiling you should make share that you have the same MPI package installed/enabled that you will
use for LAMA. 

.. code-block:: bash

   mkdir Software
   cd Software
   tar xvfz ~/Downloads parmetis-4.0.3.tar.gz
   cd parmetis-4.0.3
   export PARMETIS_ROOT=<metis_installation_directory>
   make config shared=1 cc=gcc prefix=$PARMETIS_ROOT
   make -j 8
   make install

The environment variable ``PARMETIS_ROOT`` should be used for the cmake configuration.

