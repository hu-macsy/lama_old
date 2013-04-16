Compilation and Execution of LAMA Programs on Linux Systems
===========================================================

The envrionment variable LAMA_ROOT refers the directory of your LAMA installation::

    export LAMA_ROOT=<installation/directory>

The following command compiles and links your example program simple.cpp::

    g++ -o simple simple.cpp -I${LAMA_ROOT}/include -L${LAMA_ROOT}/lib -llama 

If you have an own Boost installation, you have to add also the corresponding
include directory to the include paths::

    g++ -o simple simple.cpp -I${LAMA_ROOT}/include -I${BOOST_ROOT}/include -L${LAMA_ROOT}/lib -llama 

If this step was successful, you can run the executable::

    ./simple

Due to the dynamic linking of libraries, the executable **simple** will not contain the LAMA codes itself.
Instead, it contains a reference to the LAMA library and references will be resolved when the executable
is started. Here, it is very likely that you get the following error message::

    simple: error while loading shared libraries: liblama.so: cannot open shared object file: No such file or directory

Information about dynamically linked libraries is available by the following command::

    ldd ./simple

    linux-vdso.so.1 =>  (0x00007fff303ff000)                                                                                                    
    liblama.so => not found                                                                                                                     
    ...
    /lib64/ld-linux-x86-64.so.2 (0x00002b2c0871e000)

So the executable contains a link to the lama library but it does not know where to find this library.
There are two solutions to solve this problem.

1) Setting the library path

   You can add the lib directory to your library path. At program start, unresolved library links
   will be searched in all directories of your library path. The correspoding setting should be added
   to your bashrc file::

       export LD_LIBRARY_PATH=${LAMA_ROOT}/lib:${LD_LIBRARY_PATH}

2) Setting resolution path in the executable

   You generate a link to the LAMA lib directory within the executable. This solution is the
   preferred solution if you want to share the executable with other users::

      g++ -o simple simple.cpp -I${LAMA_ROOT}/include -L${LAMA_ROOT}/lib -llama -Wl,-rpath=${LAMA_ROOT}/lib

Now it should be possible to run the executable. Beside the output it is very likely that you get
the following warning message::

    <root> (GenLogger.cpp::275,func=configure) WARN: LAMA_LOG not set, use default configuration

The environment variable LAMA_LOG should be set with a useful value to get rid of the warning.
You can set the variable explicitly with the default value::

    export LAMA_LOG=WARN
    
For other useful environment variables see doc:`here <environmentVariables>`.

Now the output should be as follows::

    L1 norm of v = 8.8

Congratulations, your first LAMA program is running.

:download:`Download Linux Makefile <../../cpp_source/tutorial/Makefile>`

