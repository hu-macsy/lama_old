Installation
============

This description will guide you through the installation process. You find information about:

* :doc:`installation/download`: How to obtain the sources via our download portal
* :doc:`installation/requirements`: An overview of needed software (mandatory and optional)
* :doc:`installation/configuration`: CMake Configuration Details
* :doc:`installation/build`: Build (and install) step

.. Additionally you can get some tips for the installation on :doc:`windows <installation/windowsTipps>`, which is not fully tested yet.

.. toctree::
   :hidden:

   installation/download
   installation/requirements
   installation/configuration
   installation/build
   installation/windowsTipps

If your are familiar with CMake and your systems has all mandatory software prerequesites installed, all you have to do is:

.. code-block:: bash

   tar -xzvf libama-x.x.x.tar.gz
   cd libama-x.x.x
   mkdir build
   cd build
   cmake -DCMAKE_INSTALL_PREFIX=<path/to/install/dir> [options] ../scai
   make [-j <number-of-build-processes>]
