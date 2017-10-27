.. _solver-example:

Solver Example
==============

In the example directory a solver example program ``solver.exe`` can be found.
This executable can be used to test all available solvers with any input matrices.
Nearly all settings can be specified by command line arguments.

Command Line Arguments
----------------------

.. code-block:: bash

    lamaSolver.exe [options] <matrix_filename> [ rhs ] [start_solution] [ final_solution_filename ] [ options ] 

* The mandatory argument ``matrix_filename`` is the name of the file that contains the matrix.
* The optional argument ``rhs`` is either a single value or the name of the file that contains the right-hand-side.
  If not specified the right hand side is sum( matrix, dim = 2 ), i.e. each entry is the sum of the corresponding
  row. 
* The optional argument ``start_solution`` is either a single value or the name of the file that contains the intitial
  solution. If not specified the start solution is set to 0.
* If a solution filename is specified, the solution is either written into this
  file (file does not exist) or it is compared with the vector in the file (file exists).

Note: If an optional argument should be set but a previous one is not set, the previous optional argument should be ``""``.

.. code-block:: bash

    lamaSolver.exe matrix.mtx "" solution.mtx

These are the supported options:

.. code-block:: bash

         --SCAI_DISTRIBUTION=BLOCK|<dist_file_name>
         --SCAI_SOLVER=[BiCG|BiCGstab|CG|CGNE|CGNR|CGS|GMRES|InverseSolver|Jacobi|MINRES|QMR|Richardson|SimpleAMG|TFQMR]
         --SCAI_SOLVER_LOG=[noLogging|convergenceHistory|solverInformation|advancedInformation|completeInformation]
         --SCAI_MAX_ITER=<int_val>
         --SCAI_NORM=L1|L2|Max
         --SCAI_REL_TOL=<val>
         --SCAI_ABS_TOL=<val>
         --SCAI_DIV_TOL=<val>
         --SCAI_FORMAT=[CSR|ELL|JDS|DIA|COO]
         --SCAI_TYPE=[float|double|LongDouble|ComplexFloat|ComplexDouble|ComplexLongDouble]
         --SCAI_NUM_THREADS=...
         --SCAI_USE_METIS=<flag>
         --SCAI_ASYNCHRONOUS=<flag>
         or general options:
         --SCAI_COMMUNICATOR=[MPI|NO]
         --SCAI_CONTEXT=[Host|CUDA]
         --SCAI_DEVICE=[0|1|...]
         --SCAI_CUDA_USE_TEXTURE=[0|1]
         --SCAI_CUDA_USE_SHARED_MEM=[0|1]
         --SCAI_CUDA_BLOCK_SIZE=[64|128|...]

* Via ``SCAI_DISTRIBUTION`` a distribution for the matrix and the vectors can be set explicitly. 
  The value ``BLOCK`` stands for a block distribution without having an explicit mapping file.
  Otherwise the mapping is read from the corresponding file. For more information see 
  the command ``Matrix::readFromFile( matrix_filename, dist_filename)``.
* Via ``SCAI_SOLVER`` one of the available solvers can be chosen. As all solver classes register 
  themselves in a factory, they can be created dynamically at runtime.

Before the solve method is called, a summary of the whole configuration is given. Without any arguments
specified, it will print all the default settings.

.. code-block:: bash

    Solver            = CG
    Solver Logging    = convergenceHistory
    Context           = HostContext( #Threads = 8 )
    Communicator      = MPI(0:1)
    Matrix format     = CSR
    CommKind          = SYNCHRONOUS
    ValueType         = double
    #Threads/CPU      = 8
    weight            = 1
    Tolerances        = 0.00001(relative)
    Norm              = L2

Input Matrices
--------------

.. code-block:: bash

   matrix_generator.exe <filename> <dim> <stencilType> <dimX> [ <dimY> [ <dimZ> ] ]

The ``filename`` specifies the name of the output file for the matrix and the vector.

* filename = <id>.mtx -> generates matrix market format, <id>_v.mtx for vector
* filename = <id>     -> generates binary format, <id>.frm for matrix, <id>.frv for vector

``%s`` in filename is replaced with stencil values, e.g. 2D5P_100_100

The following stencil types are supported:

* dim = 1: stencilType = 3 
* dim = 2: stencilType = 5, 9
* dim = 3: stencilType = 7, 19, 27 

.. code-block:: bash

   matrix_generator matrix.mtx 3 27 50 50 50 

A Python script is provided that can be used to generate a lot of input sets for different 
sizes.

.. code-block:: bash

   python GenMatrices.py         
   ./GenMatrices.py         ! or script call

   ls data

Storage and Value Type
----------------------

You can run the solver with different value types. But be careful when reading binary data as in this
case the value type should match the type of the input data.

.. code-block:: bash

    matrix_generator.exe matrix 3 27 10 10 10           ! generates binary double data
    solver.exe matrix.frm matrix.frv                    ! okay
    solver.exe matrix.frm matrix.frv --SCAI_TYPE=float  ! not possible, binary mismatch

.. code-block:: bash

    matrix_generator matrix.mtx 3 27 10 10 10             ! generates formatted data
    solver.exe matrix.mtx matrix_v.mtx --SCAI_TYPE=float  ! okay

Solver Configuration
--------------------

The solver can be chosen by its name. The name is exactly the name that has been used
by a solver class when it has registered for the factory.

The example program `lama_info.exe` prints all registered solvers.

.. code-block:: bash

    lama_info.exe 

The :ref:`solver-logging` can be defined via ``SCAI_SOLVER_LOG``. By default, there is only
logging for the convergence history.

Every iterative solver needs a stopping criterion to define how exact the solution should be or how
many iterations the solver should run through. These criteria are connectable by using logical
connectives. For more information see :ref:`stopping-criteria`.

The example program combines the stopping criterias defined a logical ``or``, i.e. the solver 
stops if one of the criteria is fulfilled.

Parallel Execution
------------------

By default, all cores of a node are used via OpenMP parallelization in all the
matrix-vector operations used by the solver.

.. code-block:: bash

    solver.exe 3D27P_50.frm 3D27P_50.frv --SCAI_NUM_THREADS=1

On a serial machine, MPI parallelization can also be used to keep multiple cores busy.

.. code-block:: bash

    mpirun -np 2 solver.exe 3D27P_50.frm 3D27P_50.frv

As one process already uses by default all cores via OpenMP parallelization, you should
decrease the number of cores.

.. code-block:: bash

    mpirun -np 1 solver.exe 3D27P_50.frm 3D27P_50.frv --SCAI_NUM_THREADS=4
    mpirun -np 2 solver.exe 3D27P_50.frm 3D27P_50.frv --SCAI_NUM_THREADS=2
    mpirun -np 4 solver.exe 3D27P_50.frm 3D27P_50.frv --SCAI_NUM_THREADS=1

Number of threads for one process times number of processes should not be higher
than the actual number of cores available.

Of course MPI can also be used to solve a matrix on multiple nodes by using the MPI
parallelization of the matrix-vector operations.
Running MPI applications on your cluster is usually done in the following way:

.. code-block:: bash

    mpirun -machine_file machines 3D27P_50.frm 3D27P_50.frv 

Where ``machines`` is a file that contains all machines on which an MPI process is 
started.

.. code-block:: bash

    drachenfels-001
    drachenfels-002
    drachenfels-004
    drachenfels-006

When running multiple processes it can be very convenient to specifiy different configuration values
for each process:

* context, device id
* weight to indicate that load distribution should be different on the processes
* matrix storage format

.. code-block:: bash

    mpirun -machine_file machines 3D27P_50.frm 3D27P_50.frv --SCAI_WEIGHT=1,2,1,3 --SCAI_CONTEXT=Host,CUDA,Host,CUDA --SCAI_DEVICE=0,1,0,2


