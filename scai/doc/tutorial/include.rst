.. _includes:

Includes
--------

The include file lama.hpp contains some definitions how far assertions, logging and tracing statements
are compiled into your code. The definitions will be the same as used for the installation.

Usually you have to include the class definition file for each LAMA class that you are
using. As we use objects of class DenseVector and Scalar, we have to include the corresponding files.

.. TODO Explain directory structure

Data container can be found here:

.. code-block:: c++
   
   #include <scai/lama/matrices/...>        // sparse or dense matrices
   #include <scai/lama/...>                 // vectors, scalars
   #include <scai/lama/storage/...>         // sparse or dense matrix storages   

Data container related:

.. code-block:: c++
  
   #include <scai/lama/expression/...>      // Vector-Vector, Matrix-Vector or Matrix-Matrix expressions 
   #include <scai/lama/distribution/...>    // assigned to matrices and vectors for distributed use of data container
   #include <scai/lama/io/...>              // reading from / writing to container files
   #inclue  <lama/norm/...>                 // L1Norm, L2Norm, MaxNorm

Linear equation solvers:

.. code-block:: c++

   #include <scai/lama/solver/...>          // solver itself
   #include <scai/lama/solver/criteria/...> // Stopping criteria
   #include <scai/lama/solver/logger/...>   // Logging solver informations ( to file or console )

.. TODO: mpi, context and links to the API-reference
