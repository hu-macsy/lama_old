.. _Settings:

Settings
========

The class Settings provides a standardized way to query values from environment variables
and to set them via command line arguments.

Query values of environment variables
-------------------------------------

The following example shows how to query the value of the environment variable
``SCAI_SOLVER`` as a string.

.. code-block:: c++

    std::string solverName;

    if ( Settings::getEnvironment( solverName, "SCAI_SOLVER" ) )
    {
        cout << "solver = " << solverName << endl;
    }
    else
    {
        cout << "SCAI_SOLVER not defined" << endl;
    }

The following example shows how to query a boolean value from the 
environment variables ``SCAI_USE_MKL``.

.. code-block:: c++

    bool useMKL;

    if ( Settings::getEnvironment( useMKL, "SCAI_USE_MKL" ) )
    {
        cout << "useMKL = " << useMKL << endl;
    }
    else
    {
        cout << "SCAI_USE_MKL not defined or non-boolean value" << endl;
    }

A legal true value is any string starting with '1', 't', 'T', 'y', 'Y', 'j', or 'J', a legal false value
is any string starting with '0', 'f', 'F', 'n', or 'N'.

.. code-block:: c++

    SCAI_USE_MKL=0    // false
    SCAI_USE_MKL=a    // non-boolean value
    SCAI_USE_MKL=F    // false
    SCAI_USE_MKL=t    // true  

Integer values can be queried as follows:

.. code-block:: c++

    int  device;

    if ( Settings::getEnvironment( device, "SCAI_DEVICE" ) )
    {
        cout << "device = " << device << endl;
    }
    else
    {
        cout << "SCAI_DEVICE not defined" << endl;
    }

Setting environment variables via command line arguments
--------------------------------------------------------

Environment variables starting with the prefix ``SCAI_`` can also be set 
via command line arguments. This might be more convenient in certain 
situations. Therefore it is necessary to call the routine parseArgs
that takes all command line arguments starting with ``--SCAI_``

.. code-block:: c++

    int main( int argc, const char* argv[] )
    {
        common::settings::parseArgs( argc, argv );
        ...
    }

.. code-block:: c++

   SCAI variables of environment:
   ==============================
   SCAI_SOLVER=cg
   SCAI_USE_MKL=0

.. code-block:: c++

   DemoSettings --SCAI_DEVICE=7,5,1 --SCAI_SOLVER=amg 

.. code-block:: c++

   SCAI variables of environment( at program start)
   ================================================
   SCAI_SOLVER=cg
   SCAI_USE_MKL=0

   SCAI variables of environment (after parsing command line args):
   ================================================================
   SCAI_SOLVER=amg
   SCAI_USE_MKL=0
   SCAI_DEVICE=7,5,1

Especially when running parallel processes it is often desired that 
values can be set individually for each process. This is done by
using multiple values separated by a ``,``.

.. code-block:: c++

   SCAI_DEVICE=7,5,1

Getting the value of the environment variable delivers exactly one
value that is selected by the rank.

.. code-block:: c++

   Settings::setRank( 3 );  // this should be called when communicator is created
   ...
   int device = -1;
   getEnvironment( device, "SCAI_DEVICE" ) 

If the rank is greater or equal than the number n of values, the modulo operation
for the rank is taken to compute the position of the value.

.. code-block:: c++
 
   --SCAI_DEVICE=0,1,2,3,4,5     // returns device = 3 for rank = 3
   --SCAI_DEVICE=7,6,5,4         // returns device = 4 for rank = 3
   --SCAI_DEVICE=0,1             // returns device = 1 for rank = 3, pos = rank%2
   --SCAI_DEVICE=0,1,2           // returns device = 0 for rank = 3, pos = rank%3

Setting Machine Specific Environment Variables
----------------------------------------------

When running LAMA applications on heterogeneous systems it might be useful to
define environment variables that are specific for the node (processor name) and
for the rank of the process in a node. Typical examples are:

 * a variable WEIGHT might define how many load a process might take
 * a variable DOMAIN might be used to build subgroups of processors for task parallelism

For such situtations a file (e.g. ``Settings.txt``) can be used that contains such machine
dependent settings:

.. code-block:: c++

  #  comment lines and empty lines are ignored
  
  #  <machine-name> <node-rank>  <var1>=<value1> <var2>=<value2> ...

  bombelli *    WEIGHT=3.1  DOMAIN=0
  newton   0-3  WEIGHT=1.5  DOMAIN=1
  newton   4-7  WEIGHT=2.5  DOMAIN=2
  banach1  0-3  WEIGHT=4.2  DOMAIN=3
  banach2  *    WEIGHT=2.2  DOMAIN=4

The first entry in a line is the name of the node and the second entry the rank. ``*`` might be used
as a placeholder for any rank or any node. Furthermore it is possible to use ranges eg. like ``2-5``
for the rank.

Within the LAMA application the file can be read by the method ``readSettingsFile`` of the class 
Settings.

.. code-block:: c++

   const char[] settingsFileName = "Settings.txt";
   const char[] procName = ...;
   int rank = ...;
   common::Settings::readSettingsFile( settingsFileName, procName, rank );

Here is more general solution of reading a settings file whose name is specified
by the environment variable ``SCAI_SETTINGS``.

.. code-block:: c++

   std::string settingsFileName;
   const Communicator& comm = ...;    // contains processor specific settings

   if ( common::Settings::getEnvironment( settingsFileName, "SCAI_SETTINGS" ) )
   {
       common::Settings::readSettingsFile( settingsFileName.c_str(), comm.getProcessorName(), comm.getNodeRank() );
   }

The environment variables ``WEIGHT`` and ``DOMAIN`` that contain processor-specific values
might be used later to meet the corresponding requirements.

Environment Variables used in LAMA
----------------------------------

Environment variables currently used in libraries:

========================   ======  ========================================================================
Name                       Type    Description
========================   ======  ========================================================================
SCAI_LOG                   path    path to configuration file for logging
SCAI_UNSUPPORTED           string  IGNORE, WARN, or ERROR
SCAI_CONTEXT               string  specifies default context, e.g. Host, CUDA
SCAI_DEVICE                int     default device when getting a context
SCAI_THREADPOOL_SIZE       int     number of threads for asynchronous execuctions on CPU, default is 1
SCAI_ASYNCHRONOUS          int     0 (synchron), 1 (asynchron communication), or 2 (asynchron local computations) for certain routines
SCAI_CUDA_USE_CUSPARSE     bool    default true, uses cuSparse library instead of own kernels
SCAI_CUDA_USE_SHARED_MEM   bool    default true, uses cuSparse library instead of own kernels
SCAI_CUDA_USE_TEXTURE      bool    default true, uses cuSparse library instead of own kernels
SCAI_CUDA_BLOCK_SIZE       int     number of thread within one cude block, default is 128
SCAI_COMMUNICATOR          string  specifies default communicator, e.g. NO, MPI
SCAI_MPI_CUDA              bool    false, set true if MPI is CUDAaware
SCAI_USE_MKL               bool    false, use MKL library for BLAS routines
SCAI_AMG_SETUP_LIBRARY     path    library with dynamic module that that register at factory for AMG Setup
========================   ======  ========================================================================

