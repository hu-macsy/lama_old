.. _MPI:

Using MPI in LAMA
^^^^^^^^^^^^^^^^^

Message Passing Interface (MPI__) is a standardized and portable message-passing system designed by a group of researchers
from academia and industry to function on a wide variety of parallel computers. 
The standard defines the syntax and semantics of a core of library routines useful to a wide range of users writing
portable message-passing programs in Fortran 77 or the C programming language.  
LAMA uses this standard to exploit distributed-memory parallelization for operations on large matrices and vectors.

__ <https://www.mpi-forum.org/docs/docs.html>

If you want to exploit MPI in your LAMA applications, please make sure that a version of MPI (e.g. OpenMPI) is
installed on your machine, and that mpirun is in your path. By this way, the FindMPI module of cmake will recognize
your MPI installation.

LAMA itself uses the following MPI variables, which can be passed by -D VARIABLE=<path>:

 -  ``MPI_CXX_COMPILER`` for the MPI C++ compiler, e.g. mpicxx
 -  ``MPI_CXX_INCLUDE_PATH`` for the include file, e.g. /usr/local/openmpi-1.4.3/include
 -  ``MPI_CXX_LIBRARIES`` for the libraries to be linked, e.g. 
    ``/usr/local/openmpi-1.4.3/lib/libmpi_cxx.so;/usr/local/openmpi-1.4.3/lib/libmpi.so;...``

MPI is optional and LAMA can be built without it. But you will not be able to take advantage of distributed vectors and matrices.
