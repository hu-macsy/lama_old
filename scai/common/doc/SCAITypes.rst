.. _SCAITypes:

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

For every backend a list of supported arithmetic types is provided. 
E.g. the list of value types that are used in the current release:

.. code-block:: c++

   #define SCAI_ARITHMETIC_HOST     float, double, long double, ComplexFloat, ComplexDouble, ComplexLongDouble
   #define SCAI_ARITHMETIC_CUDA     float, double, ComplexFloat, ComplexDouble

For each type in this list, the templatized kernel routines will be instantiated to support
operations for vectors or matrices of the corresponding type.

The list of supported types can be configured during the CMake configuration of LAMA.
