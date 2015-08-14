:orphan:

Include
=======

.. TODO Explain directory structure

Data container can be found here:

::
   
   #include <scai/lama/matrices/...>        // sparse or dense matrices
   #include <scai/lama/...>                 // vectors, scalars
   #include <scai/lama/storage/...>         // sparse or dense matrix storages   

Data container related:

::
  
   #include <scai/lama/expression/...>      // Vector-Vector, Matrix-Vector or Matrix-Matrix expressions 
   #include <scai/lama/distribution/...>    // assigned to matrices and vectors for distributed use of data container
   #include <scai/lama/io/...>              // reading from / writing to container files
   #inclue  <lama/norm/...>            // L1Norm, L2Norm, MaxNorm

Linear equation solvers:

::

   #include <scai/lama/solver/...>          // solver itself
   #include <scai/lama/solver/criteria/...> // Stopping criteria
   #include <scai/lama/solver/logger/...>   // Logging solver informations ( to file or console )

.. TODO: mpi, context and links to the API-reference
