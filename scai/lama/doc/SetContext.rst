.. _lama_SetContext:

Setting a Context
=================

For writing algorithms you can use the generic data structures ``Scalar``, ``Vector`` and ``Matrix`` without the decision on which context the calculation should be executed, but when implementing dedicated applications you need to choose one. Otherwise the default execution context is the host and all calculation will take place on the CPU (depending on your build and environment variables: with OpenMP support).

A ``Context`` is a concept used to define where to execute a calculation and therefore where to store the data. The location of the calculation of an expression and the ``Context`` of the result is defined by the ``Context`` of the operands and the rule to copy less data to another ``Context``. Therefore in a expression with only vectors the first vector and in all other expressions the first matrix defines the execution ``Context``. All other data will be tranfered to the needed ``Context``.

The default ``Context`` is always the **Host** (CPU). Additionally, LAMA currently supports the **CUDA** context for NVIDIA GPUs. Further backends for other accelerators are planned for further releases, so you can keep up with new hardware innovations by just changing the context in your application.

You only can get an instance of a context by the factory by calling ``getContextPtr`` with a context name (and optional a device number):

.. code-block:: c++

   hmemo::ContextPtr hostCtx = hmemo::Context::getContextPtr( common::context::Host );
   hmemo::ContextPtr cudaCtx = hmemo::Context::getContextPtr( common::Context::CUDA, 0 ); // 0 defines the CUDA device used

If the requested context is not available, because it is not supported on the recent installation e.g. no CUDA is found, you receive the default context (Host). So your program using the CUDA context can be compiled and executed on another machine without CUDA running on the host without changing code. You find detailed information about ``Context``:ref:`here <scaihmemo:Context>`.

Scalar
------

A ``Scalar`` can not have context, but is always stored on the host and passed to the needed location.

Vector
------

For a ``Vector`` you can set a ``Context`` in multiple ways:


1. Set the context at creation time with the constructor (so also the initialization is done on the specific context):

.. code-block:: c++

  DenseVector<double> x( 4, 1.0, cudaCtx );

2. Set the context on a an already created vector on demand: the data associated with x will be copied to the defined context just before it is used the next time

.. code-block:: c++

  x.setContextPtr( cudaCtx ;
  
3. Set the context on a an already created vector with direct "copying": the data associated with x will be copied instantly to the defined context asynchronously (the call also returns instantly, so you can do other work while waiting for the data transfer to be done)

.. code-block:: c++

    x.prefetch( hostCtx );
    // do other usefull work
    x.wait();

Matrix
------

For a ``Matrix`` you can set a ``Context`` just on already created matrices, like 2. and 3. for a ``Vector``.

.. code-block:: c++

  CSRSparseMatrix<double> m( ... );

  // on demand "copying"
  m.setContextPtr( cudaCtx ;
 
  // direct "copying"
  m.prefetch( hostCtx );
  // do other usefull work
  m.wait();
