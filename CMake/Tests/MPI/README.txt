
OpenMPI
=======

Just make sure that mpic++, mpirun is in your PATH.
and the libraries in your LD_LIBRARY_PATH.

IntelMPI
========

Same as for OpenMPI.

Be careful about multiple MPI installations on your system. Clean your
path to have not another MPI executables in it and also take care
of the LD_LIBRARY_APTH.

Especially for Intel MPI it is difficult to detect it as it has no mpic++ which 
is preferred for the search.
