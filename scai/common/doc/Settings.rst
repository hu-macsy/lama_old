Settings
---------

 * get int, bool or string values from environment variables
 * Portable and more convenient use of ‚getenv‘ and ‚putenv‘
 * Environment variables can also be set in configuration file for logging (be careful about static initialization order)

Environment variables currently used:

========================   ======  ========================================================================
Name                       Type    Description
========================   ======  ========================================================================
SCAI_DEVICE                int     default device when getting a context
SCAI_CUDA_USE_CUSPARSE     bool    default true, uses cuSparse library instead of own kernels
SCAI_CUDA_USE_SHARED_MEM   bool    default true, uses cuSparse library instead of own kernels
SCAI_CUDA_USE_TEXTURE      bool    default true, uses cuSparse library instead of own kernels
SCAI_CUDA_BLOCK_SIZE       int     number of thread within one cude block, default is 128
SCAI_MPI_CUDA              bool    false, set true if MPI is CUDAaware
SCAI_USE_MKL               bool    false, use MKL library for BLAS routines
SCAI_AMG_SETUP_LIBRARY     path    library with dynamic module that that register at factory for AMG Setup
========================   ======  ========================================================================

Boolean variables can be defined as 0, 1 or false, true, or on, off, lower or upper case.

.. code-block:: c++

    bool useMKL;

    // using MKL for CSR might be disabled explicitly by environment variable

    if ( common::Settings::getEnvironment( useMKL, "SCAI_USE_MKL" ) )
    {
        SCAI_LOG_INFO( "useMKL = " << useMKL << ", set via SCAI_USE_MKL )
    }
    else
    {
        useMKL = true;
        SCAI_LOG_INFO( "useMKL = " << useMKL << ", is default )
    }

Instead of environment variables it is possible to use command line arguments

.. code-block:: c++

    -scai_device=0 -scai_use_mkl=off 

    int main( int argc, char** argv )
    {
        common::settings::parseArgs( argc, argv );
        ...
    }

Note that communicator must sets the node rank within its constructor.

.. code-block:: c++

  class Settings
  {
  public:

      /** Set a flag by value of its environment variable
       *
       *  @param[out]  flag is boolean variable that will be set
       *  @param[in]   envVarName is name of the environment variable
       *  @return      true if environment variable has been used to set flag
       */
      static bool getEnvironment( bool& flag, const char* envVarName );
   
      /** Set a integer by value of its environment variable */
      static bool getEnvironment( int& val, const char* envVarName );
  
      /** Set a string by value of its environment variable */
      static bool getEnvironment( std::string& val, const char* envVarName );
  
      /** Define an environment variable */

      static void putEnvironment( const char* envVarName, const char* val, bool replace = true );

      /** Specify which value has to be select for vals like <val0>:<val1> */

      static void setNodeRank( int noderank, int nodesize )
  }

setNodeRank Should be called by the communicator constructor.

Arguments for different node ranks:

.. code-block:: c++

    SCAI_DEVICE=0:2
    -scai_device=0:2

