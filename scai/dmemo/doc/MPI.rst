.. _MPI:

Using MPI in LAMA
^^^^^^^^^^^^^^^^^

The |MPI| is a standardized and portable message-passing system designed by a group of researchers
from academia and industry to function on a wide variety of parallel computers. 
The standard defines the syntax and semantics of a core of library routines useful to a wide range of users writing
portable message-passing programs in Fortran 77 or the C programming language.  
LAMA uses this standard to exploit distributed-memory parallelization for operations on large matrices and vectors.

.. |MPI| raw:: html

  <a href="https://www.mpi-forum.org/docs/docs.html" target="_blank">Message Passing Interface (MPI)</a>

If you want to exploit MPI in your LAMA applications, please make sure that a version of MPI (e.g. OpenMPI) is
installed on your machine, and that mpirun is in your path. By this way, the FindMPI module of cmake will recognize
your MPI installation.

LAMA itself uses the following MPI variables, which can be passed by -D VARIABLE=<path>:

 -  ``MPI_CXX_COMPILER`` for the MPI C++ compiler, e.g. mpicxx
 -  ``MPI_CXX_INCLUDE_PATH`` for the include file, e.g. /usr/local/openmpi-1.4.3/include
 -  ``MPI_CXX_LIBRARIES`` for the libraries to be linked, e.g. 
    ``/usr/local/openmpi-1.4.3/lib/libmpi_cxx.so;/usr/local/openmpi-1.4.3/lib/libmpi.so;...``

In contrary to many other software packages that use MPI, the build does not use any MPI compiler driver like 
``mpicxx``. This has the advantage that there will be never conflicts between the C++ compiler used by the MPI driver
and the C++ compiler specified for the CMake configuration. 

MPI is optional and LAMA can be built without it. But you will not be able to take advantage of distributed vectors and matrices.
