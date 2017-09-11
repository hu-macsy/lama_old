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

.. code-block:: c++

    export PATH=${SCAI_HOME}/bin:${PATH}

*********
Reference
*********

These are the provided tools:

=======================    ==========================================
Tool/Executable            Description
=======================    ==========================================
:ref:`LamaInfo`                
:ref:`LamaSolver`
=======================    ==========================================

.. toctree::
   :hidden:

   LamaInfo
   LamaSolver

*******
Example
*******


.. code-block:: c++


************
Dependencies
************

The tools executable are dependant of all underlying libraries:

* :ref:`SCAI LAMA <scailama:main-page_lama>`

************
Related Work
************

None.
