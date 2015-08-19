Logging
=======

Inserting log statements into code is a low-tech method for debugging it. It may also be the only way because
debuggers are not always available or applicable. This is usually the case for multithreaded applications and
distributed applications at large.

Experience indicates that logging is an important component of the development cycle. It offers several
advantages. It provides precise context about a run of the application. Once inserted into the code, the
generation of logging output requires no human intervention. Moreover, log output can be saved in persistent
medium to be studied at a later time. In addition to its use in the development cycle, a sufficiently rich
logging package can also be viewed as an auditing tool.

Logging does have its drawbacks. It can slow down an application. If too verbose, it can cause scrolling
blindness. To alleviate these concerns, our logging module is designed to be reliable, fast and extensible.
Since logging is rarely the main focus of an application, the API strives to be simple to understand and to
use.

Logging Libraries
-----------------

The following three versions of Logging libraries are available:

- `log4cpp`_
- `log4cxx`_
- `boost logging <http://boost-log.sourceforge.net/libs/log/doc/html/index.html>`_

.. _log4cpp: http://log4cpp.sourceforge.net/
.. _log4cxx: http://logging.apache.org/log4cxx/

There are certain advantages and disadvantges for each of the libraries.
There is a comparative survey `here`__ for the different facilities but this survey does not give any hints
about performance, usability, and availability. Here are some impressions for each of the libraries:

__ http://log4cpp.hora-obscura.de/index.php/LoggingLibraryForCpp

- log4cpp is very simple, easy to install, but not very performant.

- log4cxx is the most popular, very performant, highest functionality, but requires the Apache runtime system
  and utilities to be installed.
  
- Boost logging might be integrated in one of the next versions of Boost but not yet; it has only include
  files but no library itself; runtime configuration files are not supported.

SCAI Logging
------------

We decided to use an own simple logging facilitate to decrease dependencies with other software packages.
The macros make it possible to implement the logging with one arbitrary library and may be later
with an own version.

Each of the following macro starts with ``SCAI_LOG`` to indicate that it is used for logging.

Global definition, configuration and initialization of the logging system is not necessary.
It is done implicitly when the first logging statement is used.

In the different subroutines and modules it is possible to get access to a logger. The macro SCAI_LOG_DEF_LOGGER
defines a static logger variable with a certain name. By the dot notation for the name loggers can be structured
hierarchically that makes it more comfortable to configure the loggers.

::

    SCAI_LOG_DEF_LOGGER( Vector::logger, "Vector" ); 
    SCAI_LOG_DEF_LOGGER( CSRMatrix::logger, "Matrix.CSR" );
    SCAI_LOG_DEF_LOGGER( ELLMatrix::logger, "Matrix.ELL" );

Furthermore, for template classes:

::

    SCAI_LOG_DEF_TEMPLATE_LOGGER( template<typename T>, SparseMatrix<T>::logger, "Matrix.SparseMatrix" )

A logger is used in the following macros that stand for the logging statements at the different levels. The
variable logger must have been defined with one of the two previous macros:

::

    SCAI_LOG_TRACE( logger, "trace message" );
    SCAI_LOG_DEBUG( logger, "debug message" );
    SCAI_LOG_INFO ( logger, "info message" );
    SCAI_LOG_WARN ( logger, "warn message" );
    SCAI_LOG_ERROR( logger, "error message" );
    SCAI_LOG_FATAL( logger, "fatal message" );

It is possible to combine arguments like it is done for streams:

::

    SCAI_LOG_DEBUG( logger, "loop iteration " << i " of " << n );

The general idea is that the logging should appear in the source code but logging is usually disabled at
runtime especially for the lower levels (DEBUG, INFO).
The logging macros should have nearly no overhead if the logging is disabled at runtime. 
Rather good performance is achieved if there is only one integer comparison at runtime for a logging statement. 
It could be better if this comparison can even be extracted out of loops. But it could be worse if e.g. string
operations are executed for the logging messages even if logging is disabled.

There is an include file that contains the definitions for all the macros:

::

    #include <scai/logging.hpp>

Configuration logging at compile time
-------------------------------------

By default, all logging statements will be compiled so that they can be enabled at runtime. If the logging
is switched off at runtime, the logging becomes just a comparison between integers. But even this might
cause some overhead so this logging library provides the possibility to disable logging already at compile time.

Before the file "scai/logging.hpp" is included one of the following macros might be used to disable certain logging
messages in the source file already at compile time:

::

    #define SCAI_LOG_LEVEL_TRACE
    #define SCAI_LOG_LEVEL_DEBUG    // no trace logging
    #define SCAI_LOG_LEVEL_INFO     // default, no trace and debug logging
    #define SCAI_LOG_LEVEL_WARN     // no trace, debug and info logging
    #define SCAI_LOG_LEVEL_ERROR    // no trace, debug, info, warn logging
    #define SCAI_LOG_LEVEL_FATAL    // no trace, debug, info, warn, error
    #define SCAI_LOG_LEVEL_OFF      // all logging statements are disabled

Be careful when using one of these macros. It implies that the corresponding logging messages for the lower
levels will become empty and do not appear any more in the code. They cannot be enabled at runtime without
recompilation. The advantage is that it allows better optimized code as there is not even any integer
comparison in the code.

Conditional code for logging
----------------------------

In some situations it might be necessary to have some additional code that is needed to compute values for
the logging. It should be executed only if the corresponding logging level is enabled.

::

   #ifdef SCAI_LOG_INFO_ENABLED
   if ( SCAI_LOG_INFO_ON( rootLogger ) )
   {
       int sum = 0;
       for ( int k = 0; k < N; k++ )
       {
           sum += k;
       }
       SCAI_LOG_INFO( rootLogger, "main program terminates with sum = " << sum );
   }
   #endif

The macro SCAI_LOG_INFO_ON( logger ) returns true if the info level is enabled for the logger at runtme. The
guard ``SCAI_LOG_INFO_ENABLED`` might be used to ignore the code even at compile time if the LOG_LEVEL is
higher than INFO.

Use of logging for C++ classes
------------------------------

Usually, each C++ class should have its own logger that is used within the methods of the class. 
The logger becomes a static variable of the class.

::

   #include "scai/logging.hpp"
   
   class Example
   {
       ...
   protected: 
       SCAI_LOG_DECL_STATIC_LOGGER( logger )
       ...
   }

   template<typename T>
   class SparseMatrix
   {
   protected: 
       SCAI_LOG_DECL_STATIC_LOGGER( logger )
       ...
   }

As the logger is a static variable,
it is not possible to set logging levels individually for objects of the class. This seemed to be nice
but has two major problems. The first one is an efficiency reason as each construction of an object requires
not very cheap access to the logger in the logger hierarchy. The second one is that the  configuration of
loggers for individual objects is not practical as objects have no individual names.

A logger should not be declared as public so the logger can be used only inside the class. 
Derived classes should usually have their own logger, in this case
the logger should become private. The logger should be protected in situatons  where it is 
useful that the logger can also be used in derived classes.

In the implementation of the class, e.g. Example.cpp, the logger has to be defined as follows:

::

    SCAI_LOG_DEF_LOGGER( Example::logger, "Example" )
    SCAI_LOG_DEF_TEMPLATE_LOGGER( template <typename T>, SparseMatrix<T>::logger, "Matrix.SparseMatrix" )
 
Configuration of logging at runtime
-----------------------------------

Logging can be configured at runtime by setting the environment variable ``SCAI_LOG`` with a configuration file.

.. code-block:: bash

    export SCAI_LOG=config

If the variable is not set, a logging file with the name .loggingrc is searched in the home directory of the user.
If this file is also not available, the default configuration is chosen.

The configuration file should contain lines that specfy the levels of the logger.

::

    <root> = ERROR
    Matrix = INFO
    Matrix.CSRSparseMatrix = DEBUG
    Distribution = INFO
    Distribution.BlockDistribution = WARN

The default configuration for all loggers is level *WARN* if no configuration file is specified or if no
level has been specified in the configuration file. The RootLogger can be referenced by **<root>**.

For Debugging purposes it is also possible to flush the output of the logger, so all logging messages are
displayed even if the program crashes. Flushing can be activated by the config file:

::

    flush = true

The default output format of logging messages is as follows:

::

    #date, #time #name @ #thread ( #func -> #file::#line ) #level #msg

where the tokens starting with # have the following meanings:

- ``#date`` stands for the current date, e.g. 2015-07-26 (yyyy-mm-dd)
- ``#time`` stands for the time of the output, e.g. 13:21:22 (hh:mm:ss)
- ``#name`` stands for the full name of the logger
- ``#func`` stands for the function in which the logging has been called
- ``#file`` is the file contaning the logging macro
- ``#line`` is the line number in the file with the actual logging statement
- ``#level`` is the logging level (e.g. INFO or WARN)
- ``#msg`` is the output message of the logging statement

It is possible to change this default output format by a line in the config file, e.g.:

::

    format = "logger = #name, msg: #msg"

The output format cannot be redefined individually for different loggers.

Compile Flags for Logging
-------------------------

Source files containing logging statements should include the file ``scai/logging.hpp``.

::

    #include <scai/logging.hpp>

For compilation, the corresponding include directory must be set.

Furthermore, a compile flag ``SCAI_LOG_LEVEL_xxx`` must be specified to specify
which logging statement should be inserted in the code. All logging statements of level
xxx and higher are included, logging statements with a lower level be considered as deleted.
Good choices are:

- DEBUG should be chosen for DEBUG mode
- INFO should be chosen in RELEASE mode
- TRACE should be set in case of serious problems
- OFF might be used for benchmarking.

As logging does not cause much overhead when it is switched off at runtime, the DEBUG level is 
usually the first choice. The TRACE level might cause some overhead as it might be used in 
innermost loops.

Please keep in mind that setting a certain level at compile time will remove all logging statements with a
lower level and they can not be used at runtime any more.
 
Support for Logging Library with CMake
--------------------------------------

The logging library itself is built with CMake.

For using the logging library an CMake module ``Findscai_logging.cmake`` is provided.

::

    # find installation of logging library

    find_package( scai_logging REQUIRED )

    if ( SCAI_LOGGING_FOUND )
       # set the include directory containing logging.hpp
       include_directories( ${SCAI_LOGGING_INCLUDE_DIR} )
       # set compilation flag, same as -DSCAI_LOG_${SCAI_LOGGING_LEVEL} )
       add_definitions( -D${SCAI_LOGGING_FLAG} )
       ...
       target_link_libraries( <newlib> ${SCAI_LOGGING_LIBRARY} )
    endif ( SCAI_LOGGING_FOUND )

The default logging level is chosen by the built type.

::

    #  CMAKE_BUILD_TYPE == Debug   : use -DSCAI_LOG_LEVEL_DEBUG
    #  CMAKE_BUILD_TYPE == Release : use -DSCAI_LOG_LEVEL_INFO
    #
    #  For serious problems: -DSCAI_LOG_LEVEL_TRACE
    #  For benchmarks:       -DSCAI_LOG_LEVEL_OFF (or -DSCAI_LOG_LEVEL_FATAL, -DSCAI_LOG_LEVEL_ERROR)

The logging level can be set via ccmake using the CMake variable ``SCAI_LOGGING_LEVEL``.

Examples
--------

The first example shows a simple example with logging statements of all different levels.

.. literalinclude:: ../../../logging/examples/DemoLogging.cpp

At runtime it can be decided which logging will be printed:

::

    export SCAI_LOG=INFO
    examples/DemoLogging 

    2015-08-19, 16:14:11 Demo @ main ( main -> DemoLogging.cpp::44 ) INFO a message about progress in the program
    2015-08-19, 16:14:11 Demo @ main ( main -> DemoLogging.cpp::47 ) WARN a message with a warning, but execution is still possible
    2015-08-19, 16:14:11 Demo @ main ( main -> DemoLogging.cpp::48 ) ERROR a message for an error, error handling will be invoked
    2015-08-19, 16:14:11 Demo @ main ( main -> DemoLogging.cpp::49 ) FATAL a message for a fatal error, execution will stop

::

    export SCAI_LOG=ERROR
    examples/DemoLogging 

    2015-08-19, 16:14:25 Demo @ main ( main -> DemoLogging.cpp::48 ) ERROR a message for an error, error handling will be invoked
    2015-08-19, 16:14:25 Demo @ main ( main -> DemoLogging.cpp::49 ) FATAL a message for a fatal error, execution will stop

Using a config file (here with the name ``config``) is also possible, it should contain a line that sets the
level of the logger *Demo*.

::

    # config file for logging
    Demo = TRACE

::

    export SCAI_LOG=config
    2015-08-19, 16:43:25 Demo @ main ( main -> DemoLogging.cpp::44 ) INFO a message about progress in the program
    2015-08-19, 16:43:25 Demo @ main ( main -> DemoLogging.cpp::45 ) DEBUG a message useful to find bugs in the program
    2015-08-19, 16:43:25 Demo @ main ( main -> DemoLogging.cpp::46 ) TRACE a message with very detailled info, usually not compiled
    2015-08-19, 16:43:25 Demo @ main ( main -> DemoLogging.cpp::47 ) WARN a message with a warning, but execution is still possible
    2015-08-19, 16:43:25 Demo @ main ( main -> DemoLogging.cpp::48 ) ERROR a message for an error, error handling will be invoked
    2015-08-19, 16:43:25 Demo @ main ( main -> DemoLogging.cpp::49 ) FATAL a message for a fatal error, execution will stop

The next example demonstrates the use of logging with multiple threads.

.. literalinclude:: ../../../logging/examples/LogOpenMP.cpp

::

   export SCAI_LOG=INFO
   examples/LogOpenMP

   2015-08-19, 16:11:01 LogOpenMP @ OMP_Thread_4 ( main -> LogOpenMP.cpp::29 ) INFO executes iteration 8 of 20
   2015-08-19, 16:11:01 LogOpenMP @ OMP_Thread_2 ( main -> LogOpenMP.cpp::29 ) INFO executes iteration 4 of 20
   2015-08-19, 16:11:01 LogOpenMP @ OMP_Thread_4 ( main -> LogOpenMP.cpp::29 ) INFO executes iteration 9 of 20
   2015-08-19, 16:11:01 LogOpenMP @ OMP_Thread_1 ( main -> LogOpenMP.cpp::29 ) INFO executes iteration 2 of 20
   2015-08-19, 16:11:01 LogOpenMP @ OMP_Thread_3 ( main -> LogOpenMP.cpp::29 ) INFO executes iteration 6 of 20
   2015-08-19, 16:11:01 LogOpenMP @ OMP_Thread_2 ( main -> LogOpenMP.cpp::29 ) INFO executes iteration 5 of 20
   2015-08-19, 16:11:01 LogOpenMP @ OMP_Thread_5 ( main -> LogOpenMP.cpp::29 ) INFO executes iteration 10 of 20
   2015-08-19, 16:11:01 LogOpenMP @ main ( main -> LogOpenMP.cpp::29 ) INFO executes iteration 0 of 20
   2015-08-19, 16:11:01 LogOpenMP @ OMP_Thread_7 ( main -> LogOpenMP.cpp::29 ) INFO executes iteration 14 of 20
   2015-08-19, 16:11:01 LogOpenMP @ main ( main -> LogOpenMP.cpp::29 ) INFO executes iteration 1 of 20
   2015-08-19, 16:11:01 LogOpenMP @ OMP_Thread_7 ( main -> LogOpenMP.cpp::29 ) INFO executes iteration 15 of 20
   2015-08-19, 16:11:01 LogOpenMP @ main ( main -> LogOpenMP.cpp::29 ) INFO executes iteration 16 of 20
   2015-08-19, 16:11:01 LogOpenMP @ OMP_Thread_1 ( main -> LogOpenMP.cpp::29 ) INFO executes iteration 3 of 20
   2015-08-19, 16:11:01 LogOpenMP @ main ( main -> LogOpenMP.cpp::29 ) INFO executes iteration 17 of 20
   2015-08-19, 16:11:01 LogOpenMP @ OMP_Thread_1 ( main -> LogOpenMP.cpp::29 ) INFO executes iteration 18 of 20
   2015-08-19, 16:11:01 LogOpenMP @ OMP_Thread_3 ( main -> LogOpenMP.cpp::29 ) INFO executes iteration 7 of 20
   2015-08-19, 16:11:01 LogOpenMP @ OMP_Thread_1 ( main -> LogOpenMP.cpp::29 ) INFO executes iteration 19 of 20
   2015-08-19, 16:11:01 LogOpenMP @ OMP_Thread_6 ( main -> LogOpenMP.cpp::29 ) INFO executes iteration 12 of 20
   2015-08-19, 16:11:01 LogOpenMP @ OMP_Thread_5 ( main -> LogOpenMP.cpp::29 ) INFO executes iteration 11 of 20
   2015-08-19, 16:11:01 LogOpenMP @ OMP_Thread_6 ( main -> LogOpenMP.cpp::29 ) INFO executes iteration 13 of 20

Summary
-------

The SCAI logging library is a very convenient and efficient library that supports logging in C++ applications.
It has very low overhead and has quite powerful mechanisms for generating logging output.

It is used in all other SCAI projects coming with the distribution but can also be used in other applications.
A CMake module is provided that finds the include dir and the logging library and sets the compile flag
for the desired logging level.

Currently, logging is always done on standard output. It is not possible to write logging messages in different 
output streams.


