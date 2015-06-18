/**
 * @file Level.hpp
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
 * @brief Level.hpp
 * @author Thomas Brandes
 * @date 10.06.2015
 */

#pragma once

#include <ostream>
#include <string>

namespace logging
{
typedef enum
{
    TRACE, //!< even more detailed than DEBUG
    DEBUG, //!< designates fine-grained informational events
    INFO, //!< informational messages highlighting progress
    WARN, //!< for potentially harmful situations
    SERROR, //!< for errors that might still allow the application to continue
    FATAL, //!< severe errors that will presumably lead to aborts
    OFF, //!< turn logging off
    MAXLEVEL //!< end value for enumeration, also used as unknown level
} Level;

Level str2level( const std::string& value );

/** Translate level to a string.
 *
 *  \return the logging level as a string.
 */

const char* level2str( const Level level );

/** Output of level in a stream.
 */
std::ostream& operator<<( std::ostream& os, const Level& level );

}
