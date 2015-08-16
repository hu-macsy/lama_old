/**
 * @file CommonLogger.hpp
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
 * @brief Contains the CommonLogger class.
 * @author Matthias Makulla
 * @date 06.04.2011
 * @since 1.0.0
 */

#pragma once


// for dll_import
#include <scai/common/config.hpp>

// base classes
#include <scai/lama/solver/logger/Logger.hpp>

namespace scai
{

namespace lama
{

/**
 * @brief A CommonLogger which adds no prefix to its messages.
 */
class COMMON_DLL_IMPORTEXPORT CommonLogger: public Logger
{
public:
    /**
     * @brief Creates a CommonLogger with the specified properties.
     *
     * This constructor creates a CommonLogger with the specified
     * properties.
     *
     * @param[in] id             TODO[doxy] Complete Description.
     * @param[in] level          The loglevel of the logger. Messages with a
     *                           loglevel greater than the level of the logger will
     *                           be omitted. Instead of a CommonLogger with the
     *                           "noLogging" loglevel a NullLogger should be used.
     * @param[in] writeBehaviour Specifies, if the logger shall write its
     *                           output to the console and a file or to the
     *                           console/file only
     * @param[in] ignoreRank     TODO[doxy] Complete Description.
     */
    CommonLogger(
        const std::string& id,
        LogLevel::LogLevel level,
        LoggerWriteBehaviour::LoggerWriteBehaviour writeBehaviour,
        bool ignoreRank = false );

    /**
     * @brief Creates a CommonLogger with the specified properties.
     *
     * This constructor creates a CommonLogger with the specified
     * properties.
     *
     * @param[in] id             TODO[doxy] Complete Description.
     * @param[in] level          The loglevel of the logger. Messages with a
     *                           loglevel greater than the level of the logger will
     *                           be omitted. Instead of a CommonLogger with the
     *                           "noLogging" loglevel a NullLogger should be used.
     * @param[in] writeBehaviour Specifies, if the logger shall write its
     *                           output to the console and a file or to the
     *                           console/file only
     *
     * @param[in] timer          The timer which shall be used by this logger.
     * @param[in] timer          TODO[doxy] Complete Description.
     * @param[in] ignoreRank     TODO[doxy] Complete Description.
     */
    CommonLogger(
        const std::string& id,
        LogLevel::LogLevel level,
        LoggerWriteBehaviour::LoggerWriteBehaviour writeBehaviour,
        common::shared_ptr<Timer> timer,
        bool ignoreRank = false );

    /**
     * @brief Creates a CommonLogger with the specified properties.
     *
     * This constructor creates a CommonLogger with the specified
     * properties.
     *
     * @param[in] id               TODO[doxy] Complete Description.
     * @param[in] level            The loglevel of the logger. Messages with a
     *                             loglevel greater than the level of the logger will
     *                             be omitted. Instead of a CommonLogger with the
     *                             "noLogging" loglevel a NullLogger should be used.
     * @param[in] writeBehaviour  Specifies, if the logger shall write its
     *                             output to the console and a file or to the
     *                             console/file only
     * @param[in] logFileName      The filename of the file which shall be used by the
     *                             logger. May throw an exception if a logger with
     *                             a different filename has been created before.
     *
     * @param[in] timer            The timer which shall be used by this logger.
     * @param[in] ignoreRank       TODO[doxy] Complete Description.
     */
    CommonLogger(
        const std::string& id,
        LogLevel::LogLevel level,
        LoggerWriteBehaviour::LoggerWriteBehaviour writeBehaviour,
        const std::string& logFileName,
        common::shared_ptr<Timer> timer,
        bool ignoreRank = false );

    /**
     * @brief Destructor.
     */
    virtual ~CommonLogger();

    /**
     * @brief The CommonLogger adds no prefix to the message, so this method
     *        just returns an empty string.
     *
     * @return The prefix - an empty string.
     */
    std::string createPrefix();
};

} /* end namespace lama */

} /* end namespace scai */
