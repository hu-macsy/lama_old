:orphan:

Include
=======

TODO Explain directory structure

Data container can be found here:

::
   
   #include <lama/matrices/...>        // sparse or dense matrices
   #include <lama/...>                 // vectors, scalars   

Data container related:

::
  
   #include <lama/expression/...>      // Vecotr-Vector, Matrix-Vector or Matrix-Matrix expressions 
   #include <lama/distribution/...>    // assigned to matrices and vectors
   #include <lama/io/...>              // reading from / writing to container files
   #inclue  <lama/norm/...>            // L1Norm, L2Norm, MaxNorm

Linear equation solvers:

::

   #include <lama/solver/...>
   #include <lama/solver/criteria/...> // Stopping criteria
   #include <lama/solver/logger/...>   // Logging solver informations ( to file or console )

   

TODO: mpi, context and links to the API-reference