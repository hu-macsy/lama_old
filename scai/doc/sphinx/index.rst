#############
Documentation
#############

LAMA is a multi-layer framework offering four main modules (see below) for fast heterogeneous software development. 
It targets multi-core CPUs, NVIDIA GPUs and Intel Xeon Phi’s – for single-node or multi-node usage. 
LAMA’s flexible plug-in architecture allows a seamless integration of tomorrow´s CPU´s and accelerator hardware architectures, thus reducing maintenance costs.

.. figure:: _images/LAMA_Hierarchy3.png
    :width: 600px
    :align: center
    :alt: LAMA Hierarchy

The Heterogeneous Computing Development Kit provides the management of heterogeneous memory and compute kernels. Asynchronous executing is a key capability. The Math Kernel Extension gives uniform access to dense and sparse compute kernel´s on all platforms while the Distributed Extension supplies full cluster support for scalability on data parallel applications. The Linear Algebra Package enables programming using mathematical notation and prepared iterative solvers.

In the following you find a installation guide, explanation to the modules and user documentation for the containing subprojects as well as our beginners tutorial and frequently asked questions.

.. toctree::
   :titlesonly:
   :maxdepth: 1
   
   installation
   modules
   tutorial
   faq
