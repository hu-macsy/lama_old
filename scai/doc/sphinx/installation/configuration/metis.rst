Metis
^^^^^

Metis__ is a family of graph and hypergraph partitioning software. Actually we support a Metis generated (graph partitioned) Distribution, that depends on the sparse matrix structure and the number of cores. This works with *METIS_PartGraphRecursive* out of the Metis package. ParMetis comes with a full Metis installation, so you also can use it.

__ http://glaros.dtc.umn.edu/gkhome/views/metis

If your Metis installation is not in the default path define the METIS_ROOT in the cmake call with -DMETIS_ROOT=<path_to_metis> or define the environment variables METIS_INCLUDE_PATH and METIS_LIBRARY_PATH with:

.. code-block:: bash

   export METIS_INCLUDE_PATH=<path_to_metis>/include
   export METIS_LIBRARY_PATH=<path_to_metis>/lib
