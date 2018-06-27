.. _main-page_tools:

#####
Tools 
#####

***********
Description
***********

The tools module is a collection of useful exectables that also use the LAMA functionality.
These tools will be installed in the binary directory of the installation directory.

Usage
-----

It is recommended to take the bin directory in the path.

.. code-block:: bash

    export PATH=${SCAI_HOME}/bin:${PATH}

*********
Reference
*********

These are the provided tools:

============================    =======================================================
Tool/Executable                 Description
============================    =======================================================
:ref:`LamaInfo`                 Prints all registered entries of the LAMA factories
:ref:`LamaGenRandomMatrix`      Generate random matrices and write them to a file
:ref:`LamaGenStencilMatrix`     Generate stencil matrices and write them to a file
:ref:`LamaGenVector`            Generate vectors and write them to a file
:ref:`LamaMatrixConvert`        Convert a matrix file from one format to another one
:ref:`LamaVectorConvert`        Convert a vector file from one format to another one
:ref:`LamaSolver`               Framework to test different solvers
:ref:`LamaSpy`                  Generate image file to display pattern of sparse matrix
============================    =======================================================

.. toctree::
   :hidden:

   LamaInfo
   LamaGenRandomMatrix
   LamaGenStencilMatrix
   LamaGenVector
   LamaMatrixConvert
   LamaVectorConvert
   LamaSolver
   LamaSpy   

*******
Example
*******

.. code-block:: c++

    % Generate a stencil matrix
    lamaGenStencilMatrix stencil3D.frm 3 27 7 7 7
    % Generate image out.png with sparsity pattern
    lamaSpy stencil3D.frm out.png
    % Generate rhs vector
    lamaGenVector b.mtx 343 1
    % solve A * x = b with initial solution 0
    lamaSolver stencil3D.frm b.mtx 0 x.mtx

************
Dependencies
************

The tools executable are dependant of all underlying libraries:

* :ref:`SCAI LAMA <scailama:main-page_lama>`

************
Related Work
************

None.
