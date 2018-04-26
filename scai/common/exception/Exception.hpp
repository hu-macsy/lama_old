/**
 * @file Exception.hpp
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
 * @brief Definition of class Exception and macro for throwing it
 * @author Thomas Brandes
 * @date 11.06.2015
 */

#pragma once

// for dll_import
#include <scai/common/config.hpp>

// local library
#include <scai/common/NonCopyable.hpp>

// std
#include <exception>
#include <string>
#include <sstream>

namespace scai
{
/**
 * @brief The namespace common holds common stuff useful for different C++ projects
 */
namespace common
{

/**
 * @brief Exception is an exception class that contains also the call stack in its message.
 *
 * Note: default copy constructor of Exception is used when throwing an exception.
 */
class COMMON_DLL_IMPORTEXPORT Exception: public std::exception
{
public:

    /**
     * @brief The default constructor creates an Exception with no message.
     */
    Exception();

    /**
     * @brief This constructor creates an Exception with the passed message.
     *
     * @param[in] message  the message to assign to this.
     */
    Exception( const std::string& message );

    /**
     * @brief The destructor destroys this Exception.
     */
    virtual ~Exception() throw ();

    /**
     * @brief what() returns the message of this Exception.
     *
     * @return the message of this Exception.
     */
    virtual const char* what() const throw ();

    /**
     *  @brief Method that prints the current call stack in an output stream.
     *
     *  Very useful utility for identification of bugs, only supported for GNU compiler.
     */
    static void addCallStack( std::ostream& output );

protected:

    std::string mMessage;  //!< message for this exception

    /** Help routine that demangles the C++ function names of the call stack.
     *
     *  @param[in] string name of the C++ function
     *  @return demangled name as a string
     */

    static std::string demangle( const char* string );
};

}  /* end namespace common */

}  /* end namespace scai */

