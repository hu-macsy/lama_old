/**
 * @file FileLogger.hpp
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
 * @brief Contains a simple, singleton-based Logger for write only
 *        logging purposes. Designed to be used by the different Logger classes
 * @author Matthias Makulla
 * @date 06.04.2011
 */

#pragma once


// for dll_import
#include <scai/common/config.hpp>

// base classes
#include <scai/common/NonCopyable.hpp>

// std
#include <string>
#include <fstream>

namespace scai
{

namespace solver
{

/**
 * @brief Singleton-based logger for write only file-logging purposes.
 */
class COMMON_DLL_IMPORTEXPORT FileLogger: private common::NonCopyable
{
public:

    /**
     * @brief Destructor for a FileLogger object.
     */
    virtual ~FileLogger();

    /**
     * @brief Returns the only instance of the FileLogger class.
     *
     * This static method returns always the same instance of the FileLogger
     * class.
     *
     * @return The file logger
     */
    static FileLogger& getFileLogger();

    /**
     * @brief Logs a string-message.
     *
     * @param message The message which has to be logged
     */
    void logMessage( const std::string& message );

    /**
     * @brief Specifies the file name of the logfile for the logger.
     *
     * @param logFileName The name of the logfile
     * @throws Throws an exception if the name already has been set and the
     *         caller tries to set it to a different name.
     */
    void setLogFile( const std::string& logFileName );

    /**
     * @brief Closes the current logfile (if it is open)
     */
    void closeLogFile();

private:

    /**
     * @brief Filestream to the logfile.
     */
    std::fstream mFileStream;

    /**
     * @brief The name of the logfile.
     */
    std::string mFileName;

    /**
     * @brief Private constructor for the FileLogger class. In order to
     * maintain singleton characteristics this constructor has to remain
     * private
     */
    FileLogger();
};

} /* end namespace solver */

} /* end namespace scai */
