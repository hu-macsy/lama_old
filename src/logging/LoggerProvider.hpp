/**
 * @file LoggerProvider.hpp
 *
 * @license
 * Copyright (c) 2009-2013
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
 * @brief Definition of class that can create loggers.
 * @author Thomas Brandes
 * @date 02.03.2011
 * @since 1.0.0
 */
#ifndef LAMA_LOGGERPROVIDER_HPP_
#define LAMA_LOGGERPROVIDER_HPP_

// for dll import
#include <logging/config.hpp>

// others
#include <logging/AbstractLoggerCreator.hpp>

namespace log4lama
{

// creator will be used for the static loggers

AbstractLoggerCreator& theLoggerCreator();

/** Skeleton class that provides logger.
 *
 * Loggers are structured hierarchically. This class returns a reference
 * to the required logger in this hierarchy. If a logger is not available
 * it will be added before to the hierarchy.
 */
class LOG4LAMA_DLL_IMPORTEXPORT LoggerProvider
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

}

#endif // LAMA_LOGGERPROVIDER_HPP_
