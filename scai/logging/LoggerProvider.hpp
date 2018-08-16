/**
 * @file LoggerProvider.hpp
 *
 * @license
 * Copyright (c) 2009-2018
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the SCAI framework LAMA.
 *
 * LAMA is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * LAMA is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
 * more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with LAMA. If not, see <http://www.gnu.org/licenses/>.
 * @endlicense
 *
 * @brief Definition of class that can create loggers.
 * @author Thomas Brandes
 * @date 10.06.2015
 */

#pragma once

// for dll import
#include <scai/common/config.hpp>

// local library
#include <scai/logging/AbstractLoggerCreator.hpp>

namespace scai
{

namespace logging
{

// creator will be used for the static loggers

AbstractLoggerCreator& theLoggerCreator();

/** Skeleton class that provides logger.
 *
 * Loggers are structured hierarchically. This class returns a reference
 * to the required logger in this hierarchy. If a logger is not available
 * it will be added before to the hierarchy.
 */
class COMMON_DLL_IMPORTEXPORT LoggerProvider
{

public:

    virtual ~LoggerProvider();

    /**
     * @brief Returns the provider. If not yet existing, it will be created.
     * @return The provider.
     */
    static LoggerProvider& getProvider();

    /**
     * @brief Get an instance of the Logger with given id.
     * @param[in] id   The id of the logger that is needed.
     * @return Reference to the logger.
     * @throws Exception, if the id could not be resolved to a logger.
     *
     * \code
     *
     * Logger& myLogger = LoggerProvider::getProvider().getInstance("Matrix.CSR");
     * \endcode
     *
     * getInstance of a logger with the same id will return the same
     * reference. If a logger for a given id is not available yet, it
     * will be created.
     */

    class Logger& getInstance( const std::string& id ) const;

    /**
     * @brief Sets the logger creator to be used for logging.
     * @param[in] creator   The function to create the logger.
     *
     * NOTE: this routine does not help for static loggers.
     */
    void setLoggerCreator( AbstractLoggerCreator* creator );

private:

    LoggerProvider();

    mutable AbstractLoggerCreator* mLoggerCreator; //! class to create logger

    static LoggerProvider* theProvider; //!< unique instance for provider
};

} /* end namespace logging */

} /* end namespace scai */
