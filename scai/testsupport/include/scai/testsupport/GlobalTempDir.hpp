/**
 * @file include/scai/testsupport/GlobalTempDir.hpp
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
 * @brief Provides a globally shared directory for temporary files.
 * @author Andreas Longva
 * @date 21.11.2017
 */
#pragma once

#include <string>
#include <stdexcept>
#include <memory>
#include <utility>

namespace scai
{

namespace testsupport
{

/**
 * Represents a global temporary directory for use in tests.
 *
 * It is designed explicitly only for use in tests, with the following use pattern:
 *
 * 1. During test binary initialization, the temporary directory path is set.
 * 2. During execution of individual tests, the path can be acquired through getPath().
 *
 * Note that the class is NOT designed for concurrent or parallel access, and as such
 * it should only be used from a single thread. Note that the Boost test runner
 * is single-threaded, and so it is generally safe to use this class from tests.
 */
class GlobalTempDir
{
public:

    /**
     * Retrieves the temporary directory path.
     *
     * If this has not previously been set through a call to setPath(),
     * an std::logic_error exception is thrown.
     */
    inline static std::string getPath()
    {
        if ( m_tempDirPath )
        {
            return *m_tempDirPath;
        }
        else
        {
            throw std::logic_error( "Attempt to get the global temporary directory path "
                                    "before it has been initialized." );
        }
    }

    /**
     * Sets the temporary directory path.
     *
     * This can only be done exactly once throughout the lifetime of the program.
     * If a second attempt to set the path is made, an std::logic_error exception
     * is thrown.
     */
    inline static void setPath( std::string path )
    {
        if ( !m_tempDirPath )
        {
            m_tempDirPath.reset( new std::string( std::move( path ) ) );
        }
        else
        {
            throw std::logic_error( "Temporary directory path has already been set. "
                                    "It must only be set exactly once in the test binary initialization." );
        }
    }

    /**
     * Sets the temporary directory path, or if the supplied string is empty,
     * set the path to a default value instead.
     *
     * The same restrictions on usage as setPath() apply to this function.
     */
    static void setPathOrDefault( std::string path )
    {
        if ( path.empty() )
        {
            setPath( "/tmp" );
        }
        else
        {
            setPath( path );
        }
    }

private:
    // Using unique_ptr as a poor man's optional...
    static std::unique_ptr<std::string> m_tempDirPath;
};

} // namespace testsupport

} // namespace scai
