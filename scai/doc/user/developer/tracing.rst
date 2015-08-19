Tracing
=======

About Tracing
--------------

The SCAI tracing library provides the possibility to collect performance data at runtime and
to visualize it. At certain events performance counter values are read and corresponding trace
files are generated.

Source Code Instrumentation
---------------------------

A source code region in C++ code can be added as follows:

::

	SCAI_REGION( "region_name" )

Each time when the region is entered and left (end of the scope) performance counter values are read
(mainly walltime ticks) and by this way performance data for each region is collected.

The macro definition is included as follows:

::

    #include <scai/tracing.hpp>

::

    add_definitions( -DSCAI_TRACE_ON )

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

.. literalinclude:: ../../../tracing/examples/TraceCalls.cpp

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

