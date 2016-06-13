.. _blaskernel-extension:

Extension
=========

The extension of the library will be shown on an example by adding the function rot.
We will just provide the implementation for one library, the other will be handled analog. 

The following files have to be altered:

- BLASKernelTrait.hpp
- external/BLASTrait.hpp
- external/BLASWrapper.hpp
- external/BLAS_BLAS1.{hpp,cpp}

BLASKernelTrait.hpp
-------------------

The signature of rot and its id have to be added. 

.. code-block:: c++

    template<typename ValueType>
    struct rot
    {
        /**
         * @brief performs rotations of points in the plane
         *
         *        x_out = c * x_in + s * y_in
         *        y_out = c * y_in - s * x_in
         *
         * @param[in] n      number of considered elements in vectors x and y.
         * @param[in,out] x  vector with minimum ( n - 1 ) * incX + 1 elements
         * @param[in] incX   storage spacing between considered elements of x.
         * @param[in,out] y  vector with minimum ( n - 1 ) * incX + 1 elements
         * @param[in] incY   storage spacing between considered elements of y.
         * @param[in] c      scalar ( cosine of the rotation angle )
         * @param[in] s      scalar ( sine of the rotation angle )
         *
         * Note: routine has no effect if n <= 0, or incX <= 0, or incY <= 0
         */
        typedef void ( *FuncType )(
            const IndexType n,
            ValueType* x,
            const IndexType incX,
            ValueType* y,
            const IndexType incY,
            const ValueType c,
            const ValueType s );

        static const char* getId()
        {
            return "BLAS1.rot";
        }
    };

external/BLASTrait.hpp
----------------------

The class ``BLASTrait`` defines the interface to the exteranl Fortran BLAS routines.
This interface has to be extended for the four rot routines that have different names
and signatures, depending on the value type.

.. code-block:: c++

    #define CALL_DEF_ROT( type ) const scai::blaskernel::BLASTrait::BLASIndexType *, type*, const scai::blaskernel::BLASTrait::BLASIndexType *, type *, const scai::blaskernel::BLASTrait::BLASIndexType *, const type*, const type*
    
    FORTRAN_BLAS_DEF( rot, s, void, CALL_DEF_ROT( float ) );
    FORTRAN_BLAS_DEF( rot, d, void, CALL_DEF_ROT( double ) );
    FORTRAN_BLAS_DEF( rot, c, void, CALL_DEF_ROT( ComplexFloat ) );
    FORTRAN_BLAS_DEF( rot, z, void, CALL_DEF_ROT( ComplexDouble ) );
    
    #undef CALL_DEF_ROT
    
external/BLASWrapper.hpp
------------------------

The wrapper macro-defintion ``BLASWRAPPER_DEF`` has to be extended:

.. code-block:: c++

    static ValueType rot( const BLASIndexType n, ValueType *x, const BLASIndexType incX,     \ 
                          ValueType *y, const BLASIndexType incY, ValueType c, ValueType s ) \
    {                                                                                        \
        return FORTRAN_BLAS_NAME( rot, prefix1 )( &n, x, &incX, y, &incY, &c, &s );          \
    }                                                                                        \

By this way, the corresponding name of the FORTRAN BLAS routine will be used.
The macro ``FORTRAN_BLAS_NAME`` is the same that is used for the interface definitions in
the class ``BLASTrait``.

external/BLAS_BLAS1.hpp
-----------------------

The singleton class ``BLAS_BLAS1`` provides all BLAS1 wrappers and registers them in the
kernel registry.  The definition of the wrapper function for ``rot`` needs to be added in the 
class definition:

.. code-block:: c++

    /**
     * This function is the Host implementation of BLASKernelTrait::rot
     */
    template<typename ValueType>
    static void rot( const IndexType n, 
                     ValueType* x, 
                     const IndexType incX, 
                     ValueType* y,
                     const IndexType incY, 
                     const ValueType c, 
                     const ValueType s );

In contrary to the BLAS routines the wrapper routine has always the same name
as one template parameter for the ValueType is used.
        
external/BLAS_BLAS1.cpp
-----------------------

The corresponding implementation for the wrapper function ``rot`` is as follows:

.. code-block:: c++

    template<typename ValueType>
    void BLAS_BLAS1::rot( const IndexType n, 
                          ValueType* x, 
                          const IndexType incX, 
                          ValueType* y,
                          const IndexType incY, 
                          const ValueType c, 
                          const ValueType s )
    {
        SCAI_REGION( "BLAS.BLAS1.rot" )
    
        SCAI_LOG_DEBUG( logger,
           "rot<" << TypeTraits<ValueType>::id() << ">,  n = " << n << ", x = " << x << ", incX = " << incX
            << ", y = " << y << ", incY = " << incY << ", c = " << c << ", s = " << s )
    
        BLASWrapper<ValueType>::rot( static_cast<BLASTrait::BLASIndexType>( n ), 
                                     x, 
                                     static_cast<BLASTrait::BLASIndexType>( incX ),
                                     y, 
                                     static_cast<BLASTrait::BLASIndexType>( incY ), 
                                     c, 
                                     s );
    }
