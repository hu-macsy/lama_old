Compilation of Instrumented Source Code
=======================================

A flag must be set to deal correctly with the REGION macros. Tracing can be either enabled or 
disabled at compile time.

.. code-block:: bash

    add_definitions( -DSCAI_TRACE_ON )
    add_definitions( -DSCAI_TRACE_OFF )

Even if it is enabled at compile time, it must be also enabled explicitly at runtime. If
tracing has already been disabled at compile time, no tracing is possible at runtime.


CMake-Support
-------------

The tracing library itself is built with CMake.

For using the tracing library an CMake module ``Findscai_tracing.cmake`` is provided.

.. code-block:: bash

    # find installation of tracing library

    find_package( scai_tracing REQUIRED )

    if ( SCAI_TRACING_FOUND )
       # set the include directory containing tracing.hpp
       include_directories( ${SCAI_TRACING_INCLUDE_DIR} )
       add_definitions( -D${SCAI_TRACING_FLAG} )
       ...
       target_link_libraries( <newlib> ${SCAI_TRACING_LIBRARY} )
    endif ( SCAI_TRACING_FOUND )

By default, tracing is enabled. It can be disabled by the boolean CMake variable
``SCAI_TRACING``. All macros used for instrumentation will be ignored if tracing
is disabled at compile time.

The tracing library itself uses also the :ref:`SCAI logging library <scailogging:main-page_logging>`. This logging is only intended
for code development. But by setting corresponding levels of the used loggers it is possible
to see which instrumented regions are entered and left.

.. code-block:: bash

   # Tracing loggers ( section in config file of logging )
   TraceConfig = WARN               # gives warning if no runtime settings are made
   TraceRegionRecord = INFO         # logging entry at each entry and exit of a region
   TraceData = WARN
   FileTable = WARN
   RegionTable = WARN
   CallTreeTable = WARN

Using VampirTrace
-----------------

The tracing facilities might also be used as manual instrumentation for other trace file formats by calling
corresponding routines when a region is entered or left.
The Open Trace Format (OTF) was designed as a well-defined trace format with open, public domain libraries for writing and reading. 
Corresponding trace files might be visualized by tools like Vampir.

Detailed information about the Open Trace Format can be found in the |OTF| documentation.

.. |OTF| raw:: html

  <a href="https://tu-dresden.de/zih/forschung/projekte/otf" target="_blank">Open Trace Format</a>

VampirTrace itself uses this trace format for the instrumentation of parallel programs. 
It can be downloaded |VampirTrace|. The latest version 15.4.4 has been tested with LAMA.

.. |VampirTrace| raw:: html

  <a href="https://tu-dresden.de/zih/forschung/projekte/vampirtrace" target="_blank">here</a>

After download and installtion of the Vampirtrace and OTF libraries, two steps are necessary to enable this feature in LAMA:

 * For the compilation of the sources of the trace library the flag ``USE_VAMPIRTRACE`` must be enabled.
 * The corresponding VampirTrace library must be linked with SCAI tracing library.

Currently this feature is not supported via the CMake modules. Therefore it has to be added by hand in the file
``scai/tracing/CMakeLists.txt`` by replacing

.. code-block:: bash

   add_library ( ${PROJECT_NAME} ${SCAI_LIBRARY_TYPE} ${CXX_SOURCES} )

with the following entries:

.. code-block:: bash

   add_definitions ( -DUSE_VAMPIRTRACE )
   add_library ( ${PROJECT_NAME} ${SCAI_LIBRARY_TYPE} ${CXX_SOURCES} )
   target_link_libraries( ${PROJECT_NAME} ${VT_LIB_DIR}/libvt.so )
 
where ``VT_LIB_DIR`` specifies the corresponding directory where the VampirTrace library has been installed.

Hint: The Open Trace Format can also be used for MPI programs to visualize the communication between the different
processes. Therefore replace in ``scai/dmemo/CMakeLists.txt`` the corresponding line with

.. code-block:: bash

   add_library ( ${PROJECT_NAME} ${SCAI_LIBRARY_TYPE} ${CXX_SOURCES} )
   target_link_libraries( ${PROJECT_NAME} ${OTF_LIB_DIR}/libvt-mpi.so )

Attention: The development of the Open Trace Format has ended and support is not offered anymore. 
As a successor of OTF the enhanced Open Trace Format 2 (OTF2) is available in the context of the new 
Scalable Performance Measurement Infrastructure Score-P. Unfortunately this Infrastructure cannot be used with LAMA, but 
might be possible in a future LAMA release.
Nevertheless the older OTF and VampirTrace version work still with the latest Vampir release.

