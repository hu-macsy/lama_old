/**
 * @file FileLogger.hpp
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
 * @brief Contains a simple, singleton-based Logger for write only
 * logging purposes. Designed to be used by the different Logger classes
 * to make offer file logging.
 * @author Matthias Makulla
 * @date 06.04.2011
 * @since 1.0.0
 */

#pragma once


// for dll_import
#include <scai/common/config.hpp>

// base classes
#include <scai/common/NonCopyable.hpp>

#include <string>
#include <fstream>

namespace scai
{

namespace lama
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

} /* end namespace lama */

} /* end namespace scai */
