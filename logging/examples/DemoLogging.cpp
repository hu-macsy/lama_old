/**
 * @file DemoLogging.hpp
 *
 * @license
 * Copyright (c) 2009-2015
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 * @endlicense
 *
 * @brief Simple example that shows using the logging library.
 *
 * @author Thomas Brandes
 * @date 08.07.2015
 */

#include "logging/logging.hpp"

LAMA_LOG_DEF_LOGGER( myLogger, "Demo" )

int main( int argc, char** argv )
{
    // macro to give the current thread a name that appears in further logs

    LAMA_LOG_THREAD( "main" )
    
    LAMA_LOG_INFO( myLogger, "a message about progress in the program" )
    LAMA_LOG_DEBUG( myLogger, "a message useful to find bugs in the program" )
    LAMA_LOG_TRACE( myLogger, "a message with very detailled info, usually not compiled" )
    LAMA_LOG_WARN( myLogger, "a message with a warning, but execution is still possible" )
    LAMA_LOG_ERROR( myLogger, "a message for an error, error handling will be invoked" )
    LAMA_LOG_FATAL( myLogger, "a message for a fatal error, execution will stop" )
}
