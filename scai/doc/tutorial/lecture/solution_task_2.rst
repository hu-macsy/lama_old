:orphan:

Solution Task 2
---------------

Here is a solution of task 2. The code demonstrate a self-provided CG-Solver. 

.. literalinclude:: ../../../solver/examples/lecture/task2.cpp 
   :language: c++
   :lines: 35-104
   :emphasize-lines: 53-61

The emphasized code shows the self-provided algorithm of a CG-Solver.

:download:`Download complete solution Task 2 <../../../solver/examples/lecture/task2.cpp>`

**Remarks**

In the example the dense vector ``rhs`` has been initialized with a linear sequence of
values. 

.. code-block:: c++

   auto rhs = linearDenseVector<ValueType>( size, 1, 1 );

If ``ValueType`` is a complex type, the corresponding real type is used for the
representation of the residual norm. This is necessary as the comparison of complex
number is not defined.

.. code-block:: c++

    // RealType<ValueType> rnorm = norm( r );
    // RealType<ValueType> eps = 0.00001;

    auto rnorm = nomr( r );
    decltype( rnorm ) eps = 0.00001;
    ...
    for ( ... ; rnorm > eps ... )
    ....

.. csv-table::

   :header: "back to this Task", "Index", "next Task"
   :widths: 330, 340, 330

   ":doc:`task_2`", ":doc:`../lecture`", ":doc:`task_3`"
