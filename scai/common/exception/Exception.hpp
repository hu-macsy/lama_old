/**
 * @file Exception.hpp
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
    static void addCallStack( std::ostringstream& output );

protected:

    std::string mMessage;

    static std::string demangle( const char* string );
};

}  /* end namespace common */

}  /* end namespace scai */

/**
 * @brief The macro SCAI_THROWEXCEPTION throws an exception that contains
 *        source code file and line as well as call stack in its message.
 *
 * @param[in] ExceptionClass must be same or derived class from scai::common::Exception
 * @param[in] msg   message to indicate reason for the exception
 * @throws    ExceptionClass (derived from scai::common:Exception, derived from std::exception)
 */

#define SCAI_THROWEXCEPTION( ExceptionClass, msg )                             \
{                                                                              \
    std::ostringstream errorStr;                                               \
    errorStr<<"Exception in line "<<__LINE__<<" of file "<<__FILE__<<"\n";     \
    errorStr<<"    Message: "<<msg<<"\n";                                      \
    scai::common::Exception::addCallStack( errorStr );                         \
    throw ExceptionClass( errorStr.str() );                                    \
}

/** COMMON_THROWEXCEPTION just throws a simple exception */

#define COMMON_THROWEXCEPTION( msg )   \
        SCAI_THROWEXCEPTION( scai::common::Exception, msg )
