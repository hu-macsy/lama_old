First LAMA Example Program
--------------------------

The following code shows a very simple LAMA example program creating a vector of 8 entries (each with the value 1.1) and calculating and printing its L1 norm.

.. literalinclude:: ../../lama/examples/tutorial/simple.cpp 
   :language: c++
   :lines: 29-

You can :download:`download <../../lama/examples/tutorial/simple.cpp>` the source file for the next compilation and installation steps (sorry, the instructions are only for Unix systems).

Compilation
^^^^^^^^^^^

You always will need to refer to the LAMA installation directory (``SCAI_ROOT``) for the include and library path, so you should set an environment variable as abbreviation for this:

.. code-block:: bash

    export SCAI_ROOT=<installation/directory>

For the compilation of an LAMA application you have to specify
the include (``-I${SCAI_ROOT}/include``) and library (``-L${SCAI_ROOT}/lib``) path and the the highest used scai_lib (in this case: scai_lama).

If you have compiled LAMA with C++11 you have to add ``-std=c++11`` as compile flag.
Otherwise boost is used and if your Boost installation is not in the system path you also have to add also the corresponding include directory to the include paths:

So the full command for compiling and linking your example program simple.cpp looks the following in general:

.. code-block:: bash

    g++ -std=c++11 -o simple simple.cpp -I${SCAI_ROOT}/include -L${SCAI_ROOT}/lib -lscai_lama 

Note: your C++ compiler must support C++11 features that are enabled by the flag ``-std=c++11``

You may get DSO link-errors like:

.. code-block:: bash

    /usr/bin/ld: /tmp/cc4wul0P.o: undefined reference to symbol '_ZN4scai6common9PrintableD2Ev'
    /usr/bin/ld: note: '_ZN4scai6common9PrintableD2Ev' is defined in DSO [SCAI_ROOT]/lib/libscai_common.so.1.0.0 so try adding it to the linker command line
    [SCAI_ROOT]/lib/libscai_common.so.1.0.0: could not read symbols: Invalid operation
    collect2: error: ld returned 1 exit status

Then you need to define all scai libraries from the highest used scai_lib down to scai_common. When using all libs it is:

.. code-block:: bash

   -lscai_solver -lscai_lama -lscai_dmemo -lscai_sparsekernel -lscai_utilskernel -lscai_blaskernel -lscai_kregistry -lscai_hmemo -lscai_tasking -lscai_tracing -lscai_logging -lscai_common

CMake Configuration
^^^^^^^^^^^^^^^^^^^

You can set up a cmake configuration file ``CMakeLists.txt`` with the following content

.. code-block:: none

   project ( LamaExample )
   cmake_minimum_required ( VERSION 2.8 )
   find_package( SCAI REQUIRED )

   add_definitions ( ${SCAI_DEFINITIONS} )
   include_directories( ${SCAI_INCLUDE_DIRS} )
   set ( CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} ${SCAI_CXX_FLAGS} )
   add_executable ( simple simple.cpp )
   target_link_libraries( simple ${SCAI_LIBRARIES} )

The find_package command will search the package SCAI that defines the following variables:

 * ``SCAI_DEFINITIONS`` contains flags used for logging, tracing and assertions. 
 * ``SCAI_CXX_FLAGS``  contains flags like ``-std=c++11`` but also the OpenMP flag 
 * ``SCAI_INCLUDE_DIRS`` contains the above include directory, but can also contain a link to the boost include directory
 * ``SCAI_LIBRARIES``  contains the list of scai libraries

.. code-block:: bash

   cmake .

This command will probably fail as it does not find the SCAIConfig.cmake file. 

.. code-block:: bash

  Could not find a package configuration file provided by "SCAI" with any of
  the following names:

    SCAIConfig.cmake
    scai-config.cmake

.. code-block:: bash

   cmake -DSCAI_DIR=${SCAI_ROOT}

   find_package( SCAI REQUIRED PATHS ENV_${SCAI_ROOT} )


If the configuration was successful the executable simple is built just by running ``make``:

.. code-block:: bash

   make

Execution
^^^^^^^^^

If the compile step was successful, you can run the executable

.. code-block:: bash

    ./simple

Due to the dynamic linking of libraries, the executable **simple** will not contain the LAMA codes itself. 
Instead, it contains a reference to the LAMA library and references will be resolved when the executable is started. 
Here, it is very likely that you get the following error message

.. code-block:: bash

    simple: error while loading shared libraries: libscai_lama.so: cannot open shared object file: No such file or directory

Information about dynamically linked libraries is available by the following command

.. code-block:: bash

    ldd ./simple

    linux-vdso.so.1 =>  (0x00007fff303ff000)                                                                                                    
    libscai_lama.so => not found                                                                                                                     
    ...
    /lib64/ld-linux-x86-64.so.2 (0x00002b2c0871e000)

So the executable contains a link to the library libscai_lama.so but it does not know where to find this library. There are two solutions to solve this problem.

1. Setting the library path

   You can add the lib directory to your library path. At program start, unresolved library links will be searched in all directories of your library path. The correspoding setting should be added to your bashrc file:

.. code-block:: bash

       export LD_LIBRARY_PATH=${SCAI_ROOT}/lib:${LD_LIBRARY_PATH}

2. Setting resolution path in the executable

   You generate a link to the LAMA lib directory within the executable. This solution is the preferred solution if you want to share the executable with other users:

.. code-block:: bash

    g++ -std=c++11 -o simple simple.cpp -I${SCAI_ROOT}/include -L${SCAI_ROOT}/lib -lscai_lama [..] -Wl,-rpath=${SCAI_ROOT}/lib

Now it should be possible to run the executable. Beside the output it is very likely that you get the following warning message from LAMA:

.. code-block:: bash

    <root> (GenLogger.cpp::xxx,func=configure) WARN: SCAI_LOG not set, use default configuration

The environment variable SCAI_LOG should be set with a useful value to get rid of the warning.
You can set the variable explicitly with the default value

.. code-block:: bash

    export SCAI_LOG=WARN
    
You might also get a warning that the TRACE environment variable is not set.

.. code-block:: bash

    <date>, <time> TraceConfig @ thread_0 ( TraceConfig -> TraceConfig.cpp::206 ) WARN SCAI_TRACE not set, tracing is disabled. Enable by SCAI_TRACE=time|ct[:vt][:thread]

.. code-block:: bash

    export SCAI_TRACE=OFF

For other useful environment variables see :doc:`here <environmentVariables>`.

Now the output should be as follows

.. code-block:: c++

    L1 norm of v = 8.8

Congratulations, your first LAMA program is running.
