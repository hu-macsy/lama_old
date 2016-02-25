Compilation
-----------

Configuration logging at compile time
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

By default, all logging statements will be compiled so that they can be enabled at runtime. If the logging
is switched off at runtime, the logging becomes just a comparison between integers. But even this might
cause some overhead so this logging library provides the possibility to disable logging already at compile time.

Before the file "scai/logging.hpp" is included one of the following macros might be used to disable certain logging
messages in the source file already at compile time:

.. code-block:: c++

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

If more than one of the preprocessor variables is set, the one which enables the lowest level is used.

.. code-block:: c++

    #define SCAI_LOG_LEVEL_DEBUG    // debug, info, warn, error, fatal are enabled
    #define SCAI_LOG_LEVEL_OFF      // useless as lower levels are already enabled

So switching off logging by ``SCAI_LOG_LEVEL_OFF`` via define in the source code does not work
if ``-DSCAI_LOG_LEVEL_WARN`` has been set at compile flag.

Conditional code for logging
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In some situations it might be necessary to have some additional code that is needed to compute values for
the logging. It should be executed only if the corresponding logging level is enabled.

.. code-block:: c++

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


Compile Flags
^^^^^^^^^^^^^

Source files containing logging statements should include the file ``scai/logging.hpp``.

.. code-block:: c++

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


CMake-Support
^^^^^^^^^^^^^

The logging library itself is built with CMake.

For using the logging library an CMake module ``Findscai_logging.cmake`` is provided.

.. code-block:: c++

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

.. code-block:: c++

    #  CMAKE_BUILD_TYPE == Debug   : use -DSCAI_LOG_LEVEL_DEBUG
    #  CMAKE_BUILD_TYPE == Release : use -DSCAI_LOG_LEVEL_INFO
    #
    #  For serious problems: -DSCAI_LOG_LEVEL_TRACE
    #  For benchmarks:       -DSCAI_LOG_LEVEL_OFF (or -DSCAI_LOG_LEVEL_FATAL, -DSCAI_LOG_LEVEL_ERROR)

The logging level can be set via ccmake using the CMake variable ``SCAI_LOGGING_LEVEL``.