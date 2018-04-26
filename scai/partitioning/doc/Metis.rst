.. _Metis:

Metis Partitioning
^^^^^^^^^^^^^^^^^^

|METIS| is a family of graph and hypergraph partitioning software. Actually we support a Metis generated (graph partitioned) distribution, 
that depends on the sparse matrix structure and the number of cores. This works with *METIS_PartGraphRecursive* out of the Metis package. 
ParMetis comes with a full Metis installation, so you also can use it.

.. |METIS| raw:: html

  <a href="http://glaros.dtc.umn.edu/gkhome/views/metis" target="_blank">Metis</a>

If your Metis installation is not in the default path define the METIS_ROOT in the cmake call with -DMETIS_ROOT=<path_to_metis> 
or define the environment variables METIS_INCLUDE_PATH and METIS_LIBRARY_PATH as follows:

.. code-block:: bash

   export METIS_INCLUDE_PATH=<path_to_metis>/include
   export METIS_LIBRARY_PATH=<path_to_metis>/lib

Download latest stable METIS release (e.g. meits-5.1.0.tar.gz), unpack it.

.. code-block:: bash

   mkdir Software
   cd Software
   tar xvfz ~/Downloads metis-5.1.0.tar.gz
   cd metis-5.1.0
   export METIS_ROOT=<metis_installation_directory>
   make config shared=1 cc=gcc prefix=$METIS_ROOT
   make -j 8
   make install

The environment variable ``METIS_ROOT`` should be used for the cmake configuration.

