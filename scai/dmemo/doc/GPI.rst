.. _GPI:

GPI-2
^^^^^

|GPI-2| is an API for asynchronous communication. It provides a flexible, scalable and fault tolerant interface for parallel applications.

.. |GPI-2| raw:: html

  <a href="http://www.gpi-site.com/gpi2/" target="_blank">GPI-2</a>

Use of GPI in LAMA
^^^^^^^^^^^^^^^^^^

Within LAMA, GPI-2 is used for the implementation of a derived communicator class, and therefore
one possibility to run LAMA applications on distributed memory architectures.

As GPI is interoperabel with MPI, both communication libraries, i.e. both communicator
classes, MPICommunicator and GPICommunicator, can be used in one application. 

CMake Configuration for GPI
^^^^^^^^^^^^^^^^^^^^^^^^^^^

GPI-2 depends on ibverbs (Verbs library from OFED). Therefore LAMA looks for a GPI-2 and ibverbs installation to enable a build with GPI.

If your GPI-2 installation is not in the default path define the METIS_ROOT in the cmake call with -DGPI2_ROOT=<path_to_metis> or define the environment variables GPI2_INCLUDE_PATH and GPI2_LIBRARY_PATH with:

.. code-block:: bash 

   export GPI2_INCLUDE_PATH=<path_to_gpi>/include
   export GPI2_LIBRARY_PATH=<path_to_gpi>/lib

The same for ibverbs: define IBVERBS_ROOT with the cmake call or IBVERBS_INCLUDE_PATH and IBVERBS_LIBRARY_PATH in your environment.
