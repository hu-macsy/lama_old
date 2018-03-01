/**
 * @file CommonLogger.hpp
 *
 * @license
 * Copyright (c) 2009-2017
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the SCAI framework LAMA.
 *
 * LAMA is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Affero General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * LAMA is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for
 * more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with LAMA. If not, see <http://www.gnu.org/licenses/>.
 *
 * Other Usage
 * Alternatively, this file may be used in accordance with the terms and
 * conditions contained in a signed written agreement between you and
 * Fraunhofer SCAI. Please contact our distributor via info[at]scapos.com.
 * @endlicense
 *
 * @brief Contains the CommonLogger class.
 * @author Matthias Makulla
 * @date 06.04.2011
 */

#pragma once


// for dll_import
#include <scai/common/config.hpp>

// base classes
#include <scai/solver/logger/SolverLogger.hpp>

namespace scai
{

namespace solver
{

/**
 * @brief A CommonLogger which adds no prefix to its messages.
 */
class COMMON_DLL_IMPORTEXPORT CommonLogger: public SolverLogger
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
        LogLevel level,
        LoggerWriteBehaviour writeBehaviour,
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
        LogLevel level,
        LoggerWriteBehaviour writeBehaviour,
        std::shared_ptr<Timer> timer,
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
        LogLevel level,
        LoggerWriteBehaviour writeBehaviour,
        const std::string& logFileName,
        std::shared_ptr<Timer> timer,
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

} /* end namespace solver */

} /* end namespace scai */
