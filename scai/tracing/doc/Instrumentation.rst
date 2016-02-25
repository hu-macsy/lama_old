Source Code Instrumentation
===========================

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