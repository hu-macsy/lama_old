.. _kernels:

Supported BLAS and LAPACK Kernels 
=================================

The structure ``BLASKernelTrait`` contains the kernel traits for all
supported BLAS and LAPACK routines. 

Each kernel trait contains the signature of the routine and its id. 
For the BLASKernel library all kernel traits have exactly one template
argument for its ValueType. Here is the example of the kernel trait 
for the ``scal`` routine:

.. code-block:: c++

    template<typename ValueType>
    struct scal
    {
        /**
         * @brief scal replaces vector x with alpha * x.
         *
         *        x = alpha * x
         *
         * @param[in] n      number of considered elements in input vector.
         * @param[in] alpha  scalar multiplier.
         * @param[in] x      vector with minimum (n - 1) * incX + 1 elements.
         * @param[in] incX   incX storage spacing between considered elements of x.
         * @param[out] x     vector x (unchanged if n<=0 or incX <= 0 ).
         */
        typedef void ( *FuncType )(
            const IndexType n,
            const ValueType alpha,
            ValueType* x,
            const IndexType incX );

        static const char* getId()
        {
            return "BLAS1.scal";
        }
    };

Even if the usual id and the signature give already an unique key for the registry,
each id is prefixed either with ``BLAS1``, ``BLAS2``, ``BLAS3`` or ``LAPACK``.

In contrary to most other BLAS libraries, the BLAS routine has the same name for different
value types. This is very convenient for our approach as template parameters can be used
for different value types.

Corresponding to the philosophy of the kernel library, for a BLAS kernel routine all 
available implementations on the different devices can be accessed by one object.

.. code-block:: c++

    kregistry::KernelTraitContextFunction<BLASKernelTrait::scal<double> > scal;

    scal[context::Host]( ... )    // call routine on a Host device
    scal[context::CUDA]( ... )    // call routine on a CUDA device
