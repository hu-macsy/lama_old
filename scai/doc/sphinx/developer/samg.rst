:orphan:

.. _samg:

SAMG
====

This section describes how to use SAMG as an experimental feature of LAMA.

Installation
------------

The SAMG interface to the SAMG library is an experimental feature and 
not part of the public domain version. Therefore it must be downloaded separately.

.. code-block:: bash 

   svn checkout --username <developername> https://tor-2.scai.fraunhofer.de/svn/lama-amg/trunk/lama-amg lama-amg

Furthermore, you have to set a link from the LAMA source directory to this new directory.

.. code-block:: bash 

   cd <project-source>/src/lama/solver
   ln -s <lama-amg-dir>/src amg

Afterwards, the LAMA configuration will see the corresponding directory and build the interface within
the build process.

On some machines it might be necessary to use a later version of the samg library.

.. code-block:: bash 

    cd <lama-amg-dir>/src
    ldd libamg.so 
    ldd libamg_bonn.so 
    ldd libmpistubs.so 
    ldd libmpistubs_bonn.so 
    cp libamg_bonn.so libamg.so
    cp libmpistubs_bonn.so libmpistubs.so

Build and Install
-----------------

As mentioned before, due to the link to the SAMG directory, everything will be built correctly.
By the installaton, the following additional libraries will be installed:

- libamglama.so
- libamg.so
- libmpistubs.so
- libsamgplama.so
- libgpuamglama.so  (only if CUDA enabled).

In contrary to the other LAMA libraries, these libraries will not be linked against an application.
They are handled as loadable modules and will be loaded at runtime. By setting the following
environment variable, the SAMG library with its interface will be used:

.. code-block:: bash 

   export LAMA_AMG_SETUP_LIBRARY=${LAMA_ROOT}/lib/libsamgplama.so
