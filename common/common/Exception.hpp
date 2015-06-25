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
#include <common/NonCopyable.hpp>

#include <exception>
#include <string>
#include <sstream>

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

}  // namespace common

/**
 * @brief The macro COMMON_THROWEXCEPTION throws an exception that contains
 *        source code file and line as well as call stack in its message.
 *
 * @param[in] msg   message to indicate reason for the exception
 * @throws    common::Exception (derived from std::exception)
 */

#define COMMON_THROWEXCEPTION( msg )                                           \
{                                                                              \
    std::ostringstream errorStr;                                               \
    errorStr<<"Exception in line "<<__LINE__<<" of file "<<__FILE__<<"\n";     \
    errorStr<<"    Message: "<<msg<<"\n";                                      \
    common::Exception::addCallStack( errorStr );                               \
    throw common::Exception( errorStr.str() );                                 \
}

/**
 * @brief The macro COMMON_ASSERT checks a condition and throws an exception 
 *        when the condition fails.
 *
 * @param[in] cond  boolean expression that is checked
 * @param[in] msg   message to indicate reason for the exception
 * @throws    common::Exception (derived from std::exception)
 *
 * The message text will also contain the condition as string, file location
 * and the current call stack.
 */

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

/**
 * @brief The macro COMMON_ASSERT_EQUAL checks to expressions for equality
 *        and throws an exception when they are different.
 *
 * @param[in] exp1  first expression for comparison
 * @param[in] exp2  second expression for comparison
 * @param[in] msg   message to indicate reason for the exception
 * @throws    common::Exception (derived from std::exception)
 *
 * Attention: the equality operator == must be defined for the two expressions.
 *
 * The message text will also contain the expressions and their values as string, 
 * file location and the current call stack.
 */

#define COMMON_ASSERT_EQUAL( exp1, exp2, msg )                                 \
{                                                                              \
    if ( ! ( ( exp1 ) == ( exp2 ) ) )                                          \
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
