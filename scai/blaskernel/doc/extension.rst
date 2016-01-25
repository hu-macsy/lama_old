.. _blaskernel-extension:

*********
Extension
*********

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
         * @brief scal replaces vector x with alpha * x.
         *
         *        x = alpha * x
         *
         * @param[in] n      number of considered elements in input vector.
         * @param[in] x  	 vector with minimum ( n - 1 ) * incX + 1 elements
         * @param[in] incX   storage spacing between considered elements of x.
         * @param[in] y      vector with minimum ( n - 1 ) * incX + 1 elements
         * @param[in] incY   storage spacing between considered elements of y.
         * @param[in] c      scalar
         * @param[in] s      scalar
         * @param[out] x     vector x (unchanged if n<=0 or incX <= 0 ).
         * @param[out] y     vector y (unchanged if n<=0 or incX <= 0 ).
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

The interface to the function have to be definied 

.. code-block:: c++

    #define CALL_DEF_ROT( type ) const scai::blaskernel::BLASTrait::BLASIndexType *, type*, const scai::blaskernel::BLASTrait::BLASIndexType *, type *, const scai::blaskernel::BLASTrait::BLASIndexType *, const type*, const type*
    
    FORTRAN_BLAS_DEF( rot, s, void, CALL_DEF_ROT( float ) );
    FORTRAN_BLAS_DEF( rot, d, void, CALL_DEF_ROT( double ) );
    FORTRAN_BLAS_DEF( rot, c, void, CALL_DEF_ROT( ComplexFloat ) );
    FORTRAN_BLAS_DEF( rot, z, void, CALL_DEF_ROT( ComplexDouble ) );
    
    #undef CALL_DEF_ROT
    
external/BLASWrapper.hpp
------------------------

The wrapper macro-defintion has to be extended:

.. code-block:: c++

	static ValueType dot(const BLASIndexType n, ValueType *x, const BLASIndexType incX,    \ 
			ValueType *y, const BLASIndexType incY, ValueType c, ValueType s )     \
	{                                                                                      \
		return FORTRAN_BLAS_NAME( rot, prefix1 )(&n, x, &incX, y, &incY, &c, &s );     \
	}                                                       

external/BLAS_BLAS1.hpp
-----------------------

The definition of rot needs to be insert:

.. code-block:: c++

	/**
     * This function is the Host implementation of BLASKernelTrait::rot
     */
    template<typename ValueType>
    static void rot( const IndexType n, ValueType* x, const IndexType incX, ValueType* y,
    		const IndexType incY, const ValueType c, const ValueType s );
    		
external/BLAS_BLAS1.cpp
-----------------------

The corresponding implementation needs to be added:

.. code-block:: c++

	template<typename ValueType>
	void BLAS_BLAS1::rot( const IndexType n, ValueType* x, const IndexType incX, ValueType* y,
			const IndexType incY, const ValueType c, const ValueType s )
	{
	    SCAI_REGION( "BLAS.BLAS1.rot" )
	
	    SCAI_LOG_DEBUG( logger,
	    		"rot<" << TypeTraits<ValueType>::id() << ">,  n = " << n << ", x = " << x << ", incX = " << incX
	            << ", y = " << y << ", incY = " << incY << ", c = " << c << ", s = " << s )
	
	    TaskSyncToken* syncToken = TaskSyncToken::getCurrentSyncToken();
	
	    if( syncToken )
	    {
	        SCAI_LOG_WARN( logger, "asynchronous execution of BLAS1::rot not supported, do it synchronously." )
	    }
	
	    BLASWrapper<ValueType>::rot( static_cast<BLASTrait::BLASIndexType>( n ), x, static_cast<BLASTrait::BLASIndexType>( incX ),
	    		y, static_cast<BLASTrait::BLASIndexType>( incY ), c, s );
	}