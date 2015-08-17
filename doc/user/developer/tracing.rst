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

    {CMAKE_INSTALL_PREFIX}/bin/TraceCT <file.ct>

