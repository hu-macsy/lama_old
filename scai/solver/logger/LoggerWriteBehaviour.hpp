/**
 * @file LoggerWriteBehaviour.hpp
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
 * @brief Contains an enumeration which identifies whether a logger
 *        logs its messages to the console only or to a file and the console
 * @author Matthias Makulla
 * @date 06.04.2011
 * @since 1.0.0
 */
#pragma once

namespace scai
{

namespace solver
{

/**
 * @brief Contains an enumeration which identifies whether a logger
 *        logs its messages to the console only or to a file and the console
 */
namespace LoggerWriteBehaviour
{

/**
 * @brief Enumeration identifying whether a logger logs its messages to
 *        the console and a file or to the console only
 */
enum LoggerWriteBehaviour
{
    /**
     * @brief Log messages will be written to standard out only.
     */
    toConsoleOnly,

    /**
     * @brief Log messages will be written to the log file only.
     */
    toFileOnly,

    /**
     * @brief Log messages will be written to the console and the logfile.
     */
    toFileAndConsole
};

} /* end namespace LoggerWriteBehaviour */

} /* end namespace solver */

} /* end namespace scai */
