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
  files but no library itself; configuration files are not supported.

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
       for (int k = 0; k < N; k++)
       {
           sum += k;
       }
       SCAI_LOG_INFO( rootLogger, "main program terminates with sum = " << sum );
   }
	#endif

The macro SCAI_LOG_INFO_ON( logger ) returns true if the info level is enabled for the logger at runtme. The
guard LOG4_INFO_ENABLED might be used disable the code even at compile time if not needed.

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

A logger should not be declared as public. Derived classes should usually have their own logger, 
so the logger should become private. The logger should be protected in situatons  where it is 
useful that the logger can also be used in derived classes, especially if the derived class is 
a template class where no own static logger can be defined. 

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

- #date stands for the current date, e.g. 2015-07-26 (yyyy-mm-dd)
- #time stands for the time of the output, e.g. 13:21:22 (hh:mm:ss)
- #name stands for the full name of the logger
- #func stands for the function in which the logging has been called
- #file is the file contaning the logging macro
- #line is the line number in the file with the actual logging statement
- #level is the logging level (e.g. INFO or WARN)
- #msg is the output message of the logging statement

It is possible to change this default output format by a line in the config file, e.g.:

::

    format = "logger = #name, msg: #msg"

The output format cannot be redefined individually for different loggers.

Compile Flags for Logging
-------------------------

For CMake, the following variable should be set::

  SCAI_LOG_LEVEL = DEBUG ( or TRACE or INFO or OFF )

- DEBUG should be chosen for DEBUG mode
- INFO should be chosen in RELEASE mode
- TRACE should be set in case of serious problems
- OFF might be used for benchmarking.

As logging does not cause much overhead when it is switched off at runtime, the DEBUG level is 
usually the first choice. The TRACE level might cause some overhead as it might be used in 
innermost loops.

Please keep in mind that setting a certain level at compile time will remove all logging statements with a
lower level and they can not be used at runtime any more.
 
::

	#  Debug   : use -DSCAI_LOG_LEVEL_DEBUG
	#  Release : use -DSCAI_LOG_LEVEL_INFO
	#
	#  For serious problems: -DSCAI_LOG_LEVEL_TRACE
	#  For benchmarks:       -DSCAI_LOG_LEVEL_OFF (or -DSCAI_LOG_LEVEL_FATAL, -DSCAI_LOG_LEVEL_ERROR)

	ADD_DEFINITIONS( -DSCAI_LOG_LEVEL_TRACE )

Some Discussion and Further Ideas
---------------------------------

- We need some more appropriate logging levels for user output in solvers
- One idea was to set logging levels for individual objects instead of classes. This idea seemed to be nice
  but has two major problems. The first one is an efficiency reason as each construction of an object requires
  a not very cheap access to the logger in the logger hierarchy. The second one is that the  configuration of
  loggers for individual objects is not practical as objects have no individual names.
