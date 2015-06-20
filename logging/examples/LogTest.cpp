#include "logging/logging.hpp"

LAMA_LOG_DEF_LOGGER( myLogger, "LogTest" )

int main( int argc, char** argv )
{
    // macro to give the current thread a name that appears in further logs

    // LAMA_LOG_THREAD( "main" )
    
    LAMA_LOG_INFO( myLogger, "a message about progress in the program" )
    LAMA_LOG_DEBUG( myLogger, "a message useful to find bugs in the program" )
    LAMA_LOG_TRACE( myLogger, "a message with very detailled info, usually not compiled" )
    LAMA_LOG_WARN( myLogger, "a message with a warning, but execution is still possible" )
    LAMA_LOG_ERROR( myLogger, "a message for an error, error handling will be invoked" )
    LAMA_LOG_FATAL( myLogger, "a message for a fatal error, execution will stop" )
}
