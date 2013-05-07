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

Logging in LAMA
---------------

We decided to use an own simple logging facilitate to decrease dependencies with other software packages.
The macros make it possible to implement the logging of LAMA with one arbitrary library and may be later
with an own version.

Each of the following macro starts with ``LAMA_LOG`` to indicate that it is used for logging.

The first two macros are used for the global definition, configuration and initialization of the logging
system. These macro have no arguments as the usual way of configuration of the logging system is done by
reading a configuration file that itself is specified by an environment variable.

::

	LAMA_LOG_DEFINITION()
   
	int main(int argc, char **argv)
	{ // this is only done in the main file or any other initial routine
	    // it reads a configuration file for the different loggers
	    LAMA_LOG_CONFIGURE();
	    ...
	}

In the different subroutines and modules it is possible to get access to a logger. The macro LAMA_LOG_ROOTLOGGER
is used to get the root logger, with the macro LAMA_LOG_LOGGER it is possible to define a logger with a certain
name. By the dot notation for the name loggers can be structured hierarchically that makes it more
comfortable to configure the loggers.

::

	LAMA_LOG_ROOTLOGGER( rootLogger );
	LAMA_LOG_LOGGER( logger1, "Vector" ); 
	LAMA_LOG_LOGGER( logger2, "Matrix.CSR" );

A logger is used in the following macros that stand for the logging statements at the different levels. The
variable logger must have been defined with one of the two previous macros:

::

	LAMA_LOG_DEBUG(logger, "debug message");
	LAMA_LOG_INFO (logger, "info message");
	LAMA_LOG_WARN (logger, "warn message");
	LAMA_LOG_ERROR(logger, "error message");
	LAMA_LOG_FATAL(logger, "fatal message");

It is possible to combine arguments like it is done for streams:

::

	LAMA_LOG_DEBUG(logger, "loop iteration " << i " of " << n);

The general idea is that the logging should appear in the source code but logging is usually disabled at
runtime especially for the lower levels (DEBUG, INFO).
The logging macros should have nearly no overhead if the logging is disabled at runtime. 
Rather good performance is achieved if there is only one integer comparison at runtime for a logging statement. 
It could be better if this comparison can even be extracted out of loops. But it could be worse if e.g. string
operations are executed for the logging messages even if logging is disabled.

There is an include file that contains the definitions for all the macros:

::

	#include "logging/logging.hpp"

Disabling logging at compile time
---------------------------------

By default, all logging statements will be compiled so that they can be enabled at runtime. If the logging
is switched off at runtime, the logging becomes just a comparison between integers. But even this might
cause some overhead so log4espp provides the possibility to disable logging already at compile time.

Before the file "log4espp.h" is included one of the following macros might be used to disable certain logging
messages in the source file already at compile time:

::

	#define LAMA_LOG_LEVEL_TRACE
	#define LAMA_LOG_LEVEL_DEBUG    // no trace logging
	#define LAMA_LOG_LEVEL_INFO     // default, no trace and debug logging
	#define LAMA_LOG_LEVEL_WARN     // no trace, debug and info logging
	#define LAMA_LOG_LEVEL_ERROR    // no trace, debug, info, warn logging
	#define LAMA_LOG_LEVEL_FATAL    // no trace, debug, info, warn, error
	#define LAMA_LOG_LEVEL_OFF      // all logging statements are disabled

Be careful when using one of these macros. It implies that the corresponding logging messages for the lower
levels will become empty and do not appear any more in the code. They cannot be enabled at runtime without
recompilation. The advantage is that it allows better optimized code as there is not even any integer
comparison in the code.

Conditional code for logging
----------------------------

In some situations it might be necessary to have some additional code that is needed to compute values for
the logging. It should be executed only if the corresponding logging level is enabled.

::

   #ifdef LAMA_LOG_INFO_ENABLED
   if ( LAMA_LOG_INFO_ON( rootLogger ) )
   {
       int sum = 0;
       for (int k = 0; k < N; k++)
       {
           sum += k;
       }
       LAMA_LOG_INFO( rootLogger, "main program terminates with sum = " << sum );
   }
	#endif

The macro LAMA_LOG_INFO_ON( logger ) returns true if the info level is enabled for the logger at runtme. The
guard LOG4_INFO_ENABLED might be used disable the code even at compile time if not needed.

Use of logging for C++ classes
------------------------------

Usually, each C++ class should have its own logger that is used within the methods of the class. 
The logger becomes a static variable of the class.

::

   #include "logging/logging.hpp"
   
   class Example
   {
       ...
   protected: 
       LAMA_LOG_DECL_STATIC_LOGGER(logger);
       ...
   }

A logger should not be declared as public. Derived classes should usually have their own logger, 
so the logger should become private. The logger should be protected in situatons  where it is 
useful that the logger can also be used in derived classes, especially if the derived class is 
a template class where no own static logger can be defined. 

In the implementation of the class, e.g. Example.cpp, the logger has to be defined as follows:

::

	LAMA_LOG_DEF_LOGGER(Example::logger, "Example");
 
Configuration of Logging with the default logger
------------------------------------------------

Logging can be configured at runtime by setting the environment variable ``LAMA_LOG`` with a configuration file.

.. code-block:: bash

	export LAMA_LOG=config

The file config contains lines that specfy the levels of the logger.

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

Compile Flags for Logging
-------------------------

For CMake, the following variable should be set::

  LAMA_LOG_LEVEL = DEBUG ( or TRACE or INFO or OFF )

- DEBUG should be chosen for DEBUG mode
- INFO should be chosen in RELEASE mode
- TRACE should be set in case of serious problems
- OFF should be used only for benchmarking.

As logging does not cause much overhead when it is switched off at runtime, the DEBUG level is 
usually the first choice. The TRACE level might cause some overhead as it might be used in 
innermost loops.

Please keep in mind that setting a certain level at compile time will remove all logging statements with a
lower level and they can not be used at runtime any more.
 
::

	#  Debug   : use -DLAMA_LOG_LEVEL_DEBUG
	#  Release : use -DLAMA_LOG_LEVEL_INFO
	#
	#  For serious problems: -DLAMA_LOG_LEVEL_TRACE
	#  For benchmarks:       -DLAMA_LOG_LEVEL_OFF (or -DLAMA_LOG_LEVEL_FATAL, -DLAMA_LOG_LEVEL_ERROR)
	ADD_DEFINITIONS( -DLOG_LEVEL_TRACE )

Some Discussion and Further Ideas
---------------------------------

- We need some more appropriate logging levels for user output in solvers
- One idea was to set logging levels for individual objects instead of classes. This idea seemed to be nice
  but has two major problems. The first one is an efficiency reason as each construction of an object requires
  a not very cheap access to the logger in the logger hierarchy. The second one is that the  configuration of
  loggers for individual objects is not practical as objects have no individual names.
