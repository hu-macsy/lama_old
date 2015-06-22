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
#include <common/config.hpp>

#include <exception>
#include <string>
#include <sstream>

/**
 * @brief The namespace common holds common stuff useful for different C++ projects
 */
namespace common
{

/**
 * @brief The class Exception represents a general exception in common library.
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

}  // namespace common

#define COMMON_THROWEXCEPTION( msg )                                           \
{                                                                              \
    std::ostringstream errorStr;                                               \
    errorStr<<"Exception in line "<<__LINE__<<" of file "<<__FILE__<<"\n";     \
    errorStr<<"    Message: "<<msg<<"\n";                                      \
    common::Exception::addCallStack( errorStr );                               \
    throw common::Exception( errorStr.str() );                                 \
}

#define COMMON_ASSERT( cond, msg )                                             \
{                                                                              \
    if (!(cond))                                                               \
    {                                                                          \
        std::ostringstream errorStr;                                           \
        errorStr << "Assertion failed in line " << __LINE__;                   \
        errorStr << " of file " << __FILE__ << "\n";                           \
        errorStr << "    Condition: " << #cond << "\n";                        \
        errorStr << "    Message: " << msg << "\n";                            \
        common::Exception::addCallStack( errorStr );                           \
        throw common::Exception( errorStr.str() );                             \
    }                                                                          \
}

#define COMMON_ASSERT_EQUAL( exp1, exp2, msg )                                 \
{                                                                              \
    if ( ( exp1 ) != ( exp2 ) )                                                \
    {                                                                          \
        std::ostringstream errorStr;                                           \
        errorStr << "Assert equal failed in line " << __LINE__;                \
        errorStr << " of file " << __FILE__ << "\n";                           \
        errorStr << "    Message: " << msg << "\n";                            \
        errorStr << "    exp_1: " << #exp1 " = " << exp1 << "\n";              \
        errorStr << "    exp_2: " << #exp2 " = " << exp2 << "\n";              \
        common::Exception::addCallStack( errorStr );                           \
        throw common::Exception( errorStr.str() );                             \
    }                                                                          \
}
