SCAI Types
==========

The file SCAITypes.hpp includes all relevant type definitions needed 
for the instantiation of LAMA vectors and matrices.

The IndexType specifies the data type used for indexing of arrays.

.. code-block:: c++

  typedef int IndexType;      // if used arrays have less than 2^31 entries
  typedef int64_t IndexType;  // if used arrays might have more than 2^31 entries

The data type int is the preferred one as it is the same type used in external libraries
like MKL or cuSPARSE. Using ''unsigned int'' is possible but might result in a lot of
warnings at compile time.

For every backend a list of value types which are used to instantiate template is provided. 
E.g. the list of value types for the Host-Backend. 

.. literalinclude:: ../SCAITypes.hpp 
   :language: c++
   :lines: 93-109

The reason for this syntax is the use of BOOST_PP_MACRO for instantiation of 
classes and methods.

.. code-block:: c++

   // instantiate CSRStorage<ValueType> for supported arithmetic types

   #define LAMA_CSR_STORAGE_INSTANTIATE(z, I, _)                                    \  
       template class COMMON_DLL_IMPORTEXPORT CSRStorage<ARITHMETIC_HOST_TYPE_##I> ; 

   BOOST_PP_REPEAT( ARITHMETIC_HOST_TYPE_CNT, LAMA_CSR_STORAGE_INSTANTIATE, _ )

   #undef LAMA_CSR_STORAGE_INSTANTIATE                              
