/*
 * @file scai/logging/test/simpleLogging.cpp
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
 * @brief simple executable that uses each log level of SCAI logging
 * @author Jan Ecker
 * @date 03.09.2015
 * @since 2.0.0
*/

#include <scai/logging.hpp>

SCAI_LOG_DEF_LOGGER( logger, "Test" )

int main()
{
    SCAI_LOG_TRACE( logger, "trace message" )
    SCAI_LOG_DEBUG( logger, "debug message" )
    SCAI_LOG_INFO( logger, "info message" )
    SCAI_LOG_WARN( logger, "warn message" )
    SCAI_LOG_ERROR( logger, "error message" )
    SCAI_LOG_FATAL( logger, "fatal message" )

} 
