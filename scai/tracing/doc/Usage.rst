Usage
=====

Collection at runtime
---------------------

The environment variable SCAI_TRACE specifies which trace data is collected at runtime.

.. code-block:: bash

    SCAI_TRACE=time
    SCAI_TRACE=ct
    SCAI_TRACE=vt
    SCAI_TRACE=time:thread:ct

If tracing is disabled the overhead for each region is very low (just comparison with a global variable).

Timing of Regions
-----------------

.. code-block:: bash 

    export SCAI_TRACE=time
    export SCAI_TRACE=time:thread

The generated output file is called <executable>.time and contains the inclusive and exclusive costs for each region. The inclusive costs are costs for the total region, for the exclusive costs the inclusive costs of all called regions within
the region are subtracted.

.. code-block:: bash

    export SCAI_TRACE=time:PREFIX=myTest

Instead of the name of the executable the PREFIX value is used in the filename, i.e. myTest.time here.

Calltree
--------

.. code-block:: bash

    export SCAI_TRACE=ct
    export SCAI_TRACE=ct:thread  

    ${CMAKE_INSTALL_PREFIX}/bin/TraceCT <file.ct>

VampirTrace
-----------

.. code-block:: bash

    export SCAI_TRACE=vt
    ... ! run application, eg. mpirun -np 4 myTest.exe
    vampir <file>.otf

