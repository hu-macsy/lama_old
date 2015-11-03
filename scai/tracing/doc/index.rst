.. _main-page:

Tracing
=======

About Tracing
--------------

The SCAI tracing library provides the possibility to collect performance data at runtime and
to visualize it. At certain events performance counter values are read and corresponding trace
files are generated.

Source Code Instrumentation
---------------------------

A user defined region in the code allows to get performance data for an exactly defined
part of the code.
Each time when the region is entered and left (end of the scope) performance counter values are read
(mainly walltime ticks) and by this way performance data for each region is collected.

A sequence of statements can be instrumented as a region as follows:

.. code-block:: c++

    #include <scai/tracing.hpp>
    ...
    SCAI_REGION_START( "region_name" )
    <sequence_of_statements>
    SCAI_REGION_END( "region_name" )
    ...

The block of statements surrounded by these two macros should not have other exit points (e.g. return,
goto, throw). The exit points might also be instrumented by  ``SCAI_REGION_END`` but this is not
very convenient. Mismatching calls of START and END result in serious errors for collecting performance
data.

Therefore in C++ a source code region should be defined as follows:

.. code-block:: c++

    #include <scai/tracing.hpp>
    ...
    {
        SCAI_REGION( "region_name" )
        ...
    }

The macro ``SCAI_REGION`` creates an object whose constructor calls the START routine and
the destructor will call the corresponding END routine. The destructor is called automatically 
when the corresponding scope is left.

Name of Regions
---------------

Regions can be structured hierarchically by using the dot nation.

::

    SCAI_REGION( "CUDA.CSRUtils.hasDiagonalProperty" )
    SCAI_REGION( "CUDA.CSRUtils.CSR2CSC" )
    SCAI_REGION( "CUDA.CSR.matrixMultiplySizes" )
    SCAI_REGION( "CUDA.CSR.matrixAdd" )
    SCAI_REGION( "Mat.Dense.invertCyclic" )
    SCAI_REGION( "Mat.Sp.syncLocal" )
    SCAI_REGION( "Vec.Times.Mat.others" )
    SCAI_REGION( "Communicator.MPI.scatter" )

For the visualization of performance data, some tools take the first name to group regions.

Even if not supported yet, this notation might be used to explicitly configure for which regions
the collection of performance data might be enabled or disabled.

Compilation of Instrumented Source Code
---------------------------------------

A flag must be set to deal correctly with the REGION macros. Tracing can be either enabled or 
disabled at compile time.

::

    add_definitions( -DSCAI_TRACE_ON )
    add_definitions( -DSCAI_TRACE_OFF )

Even if it is enabled at compile time, it must be also enabled explicitly at compile time.

Collection at runtime
----------------------

The environment variable SCAI_TRACE specifies which trace data is collected at runtime.

::

    SCAI_TRACE=time
    SCAI_TRACE=ct
    SCAI_TRACE=time:thread:ct

If tracing is disabled the overhead for each region is very low (just comparison with a global variable).

Timing of Regions
-----------------

:: 

    export SCAI_TRACE=time
    export SCAI_TRACE=time:thread

The generated output file is called <executable>.time and contains the inclusive and exclusive costs for each region. The inclusive costs are costs for the total region, for the exclusive costs the inclusive costs of all called regions within
the region are subtracted.

:: 

    export SCAI_TRACE=time:PREFIX=myTest

Instead of the name of the executable the PREFIX value is used in the filename, i.e. myTest.time here.

Calltree
--------

:: 

    export SCAI_TRACE=ct
    export SCAI_TRACE=ct:thread  

    ${CMAKE_INSTALL_PREFIX}/bin/TraceCT <file.ct>

Example
-------

The following example shows a C++ program where the main program and the two subroutines have been
instrumented.

.. literalinclude:: ../examples/TraceCalls.cpp

::

   export SCAI_TRACE=ct
   ./TraceCalls
   $HOME/local/scai/bin/TraceCT TraceCalls.ct

In the callgraph each node stands for a region and an edge from region_1 to region_2 
indicates that region_2 has been called at least once from region_1.
In the first image, the nodes and edges are labeled with the time spent in the different regions.
The values are relative, i.e. range from 0 to 100%. The left value in a node is the inclusive
time and the right value is the exclusive time. This is for leaf nodes always the same. For the main
program the inclusive time is always 100%, but in the routine only 77.21% of the time has been 
spent. The other time has been spent in the routines A and B. All exclusive times sum always up
to 100% ( 77.21 + 13.82 + 13.97 = 100 ).

.. figure:: /_images/CallTreeEvent.png
    :width: 500px
    :align: center
    :alt: Event Costs

The property window allows to set many options for the output image. Now we select
MetricCosts, where ticks are scaled to milliseconds.

.. figure:: /_images/CallTreeProps.png
    :width: 500px
    :align: center
    :alt: Event Costs

After pushing the Apply buttion a new image is generated.

.. figure:: /_images/CallTreeTime.png
    :width: 500px
    :align: center
    :alt: Event Costs

Runtime calls show how often a subroutine has been called.

.. figure:: /_images/CallTreeCalls.png
    :width: 500px
    :align: center
    :alt: Runtime Calls

Support for Tracing Library with CMake
--------------------------------------

The tracing library itself is built with CMake.

For using the tracing library an CMake module ``Findscai_tracing.cmake`` is provided.

::

    # find installation of tracing library

    find_package( scai_tracing REQUIRED )

    if ( SCAI_TRACING_FOUND )
       # set the include directory containing tracing.hpp
       include_directories( ${SCAI_TRACING_INCLUDE_DIR} )
       add_definitions( -D${SCAI_TRACING_FLAG} )
       ...
       target_link_libraries( <newlib> ${SCAI_TRACING_LIBRARY} )
    endif ( SCAI_TRACING_FOUND )

By default, tracing is enabled. It can be disabled by the boolean CMake variable
``SCAI_TRACING``. All macros used for instrumentation will be ignored if tracing
is disabled at compile time.

The tracing library itself uses also the :ref:`SCAI logging library <scailogging:main-page>`. This logging is only intended
for code development. But by setting corresponding levels of the used loggers it is possible
to see which instrumented regions are entered and left.

::

   # Tracing loggers ( section in config file of logging )
   TraceConfig = WARN               # gives warning if no runtime settings are made
   TraceRegionRecord = INFO         # logging entry at each entry and exit of a region
   TraceData = WARN
   FileTable = WARN
   RegionTable = WARN
   CallTreeTable = WARN

