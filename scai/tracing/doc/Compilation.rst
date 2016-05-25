Compilation of Instrumented Source Code
=======================================

A flag must be set to deal correctly with the REGION macros. Tracing can be either enabled or 
disabled at compile time.

.. code-block:: bash

    add_definitions( -DSCAI_TRACE_ON )
    add_definitions( -DSCAI_TRACE_OFF )

Even if it is enabled at compile time, it must be also enabled explicitly at compile time.


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

