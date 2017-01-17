.. _GPI:

GPI-2
^^^^^

|GPI-2| is an API for asynchronous communication. It provides a flexible, scalable and fault tolerant interface for parallel applications.

.. |GPI-2| raw:: html

  <a href="http://www.gpi-site.com/gpi2/" target="_blank">GPI-2</a>

Use of GPI-2 in LAMA
^^^^^^^^^^^^^^^^^^^^

Within LAMA, GPI-2 is used for the implementation of a derived communicator class, and therefore
one possibility to run LAMA applications on distributed memory architectures.

Currently, version GPI-2-1.3.0 or higher is supported.

As proposed in the GASPI standard, GPI-2 is interoperabel with MPI.
Therefore, both communication libraries, i.e. both communicator
classes, MPICommunicator and GPICommunicator, can be used in one application. 

CMake Configuration for GPI
^^^^^^^^^^^^^^^^^^^^^^^^^^^

The GPI-2 implementaion based on Infiniband/RoCE requires ibverbs (Verbs library from OFED). 
Therefore cmake searches for a GPI-2 and ibverbs installation to enable a build with GPI.
GPI-2 can also be configured in such a way that it uses TCP sockets and the ethernet
without the need of Infiniband and ibverbs.

If your GPI-2 installation is not in the default path, define the ``GPI2_ROOT`` in the cmake call with ``-DGPI2_ROOT=<path_to_gpi2>``.
or set the corresponding environment variable:

.. code-block:: bash 

   export GPI2_ROOT=<path_to_gpi>
   cmake ...

Setting the root variable has the same effect as setting the both variables 
``GPI2_INCLUDE_PATH`` and ``GPI2_LIBRARY_PATH`` as follows:

.. code-block:: bash 

   export GPI2_INCLUDE_PATH=$GPI2_ROOT/include
   export GPI2_LIBRARY_PATH=$GPI2_ROOT/lib

define separately the environment variables ``GPI2_INCLUDE_PATH`` and ``GPI2_LIBRARY_PATH`` with:

The same for ibverbs: define IBVERBS_ROOT with the cmake call or IBVERBS_INCLUDE_PATH and IBVERBS_LIBRARY_PATH in your environment.

Restrictions for GPI
^^^^^^^^^^^^^^^^^^^^

The current implementation for the GPI communicator cannot support the value type ``ComplexLongDouble``.
The reason is that the memory management of the GASPI segments allocates data on 16-Byte boundaries and so
there are serious alignment problems for 32-Byte data values.

.. code-block:: bash 

    SCAI_HOST_TYPES_LIST=float;double;long;ComplexFlag;ComplexDouble

