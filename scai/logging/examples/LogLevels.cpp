/**
 * @file LogLevels.cpp
 * @brief Simple example that shows using the logging library.
 */

#include <scai/logging.hpp>

SCAI_LOG_DEF_LOGGER( myLogger, "Demo" )

int main( int, char** )
{
    // macro to give the current thread a name that appears in further logs

    SCAI_LOG_THREAD( "main" )
    
    SCAI_LOG_INFO( myLogger, "a message about progress in the program" )
    SCAI_LOG_DEBUG( myLogger, "a message useful to find bugs in the program" )
    SCAI_LOG_TRACE( myLogger, "a message with very detailled info, usually not compiled" )
    SCAI_LOG_WARN( myLogger, "a message with a warning, but execution is still possible" )
    SCAI_LOG_ERROR( myLogger, "a message for an error, error handling will be invoked" )
    SCAI_LOG_FATAL( myLogger, "a message for a fatal error, execution will stop" )
}
