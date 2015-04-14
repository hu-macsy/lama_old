/**
 * @file LAMAAssert.hpp
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
 * @brief Contains error checking macros like ASSERT and LAMACALL.
 * @author Jiri Kraus
 * @date 02.03.2011
 * @since 1.0.0
 */

#ifndef LAMA_ASSERT_HPP_
#define LAMA_ASSERT_HPP_

#include <lama/exception/Exception.hpp>

#include <sstream>
#include <cstdio>
#include <iostream>

//Little help for the Eclipse Parser
#ifdef __CDT_PARSER__
#define LAMA_ASSERT_LEVEL_DEBUG
#endif

#ifndef LAMA_ASSERT_LEVEL_OFF
#    ifndef LAMA_CHECK_ASSERTS
/**
 * @brief The macro LAMACHECKASSERTS is used to control the checking of
 *        asserts.
 *
 * The macro LAMACHECKASSERTS is used to control the checking of asserts.
 * LAMACHECKASSERTS will be automatically defined if NDEBUG (for debug
 * builds) is not defined.
 * If LAMACHECKASSERTS is not defined the assertions will not be checked.
 */
#        define LAMA_CHECK_ASSERTS
#    endif
#endif // NDEBUG
/** Dummy function to avoid compiler warnings. */

template<typename ValueType>
static inline void unused( ValueType const )
{
    // do nothing
}

template<typename ValueType1,typename ValueType2>
static inline void unused( const ValueType1&, const ValueType2& )
{
    // do nothing
}

#ifndef LAMA_CHECK_ASSERTS

/**
 * @brief LAMA_ASSERT(exp, msg) is a more sophisticated assert macro.
 *
 * LAMA_ASSERT(exp, msg) is a more sophisticated assert macro. If LAMACHECKASSERTS is
 * defined and the expression exp does not evaluate to true an lama::Exception
 * is thrown. The macro will generate a message for the exception with this
 * content
 * "Assertion Failed in line <lineno> of file <filename>\n\tMessage: msg\n"
 * As the message for the exception is build using a stringstream msg can be
 * everything that can be appended to an output stream, e.g.
 * " The dimension (="<<v.size()<<") needs to be at least "<<minSize
 * If LAMACHECKASSERTS is not defined LAMA_ASSERT(exp,msg) will be replaced by nothing.
 *
 * @param exp   the expression to test. [IN]
 * @param msg   the message to assign to the exception. [IN]
 * @throws      lama::Exception
 */
#define LAMA_ASSERT(exp, msg) unused( exp );

#define LAMA_ASSERT_EQUAL(exp1, exp2) unused( exp1, exp2 );

/**
 * @brief LAMA_CHECK_ERROR(msg) is a convenience macro to check the last error
 *        that has been occurred in the lama library.
 *
 * CHECKLAMAERROR(msg) is a convenience macro to check the last error that has
 * that has been occurred in the lama library. If the last error is not
 * LAMA_STATUS_SUCCESS an lama::Exception is thrown. It will be build like
 * the exception thrown by the macro LAMACALL.
 *
 * If LAMA_CHECK_ASSERTS is not defined LAMA_CHECK_ERROR(msg) will not check the last
 * error.
 *
 * @param msg   the message to assign to the exception. [IN]
 * @throws      lama::Exception
 */
#define LAMA_CHECK_ERROR(msg)

#else // LAMACHECKASSERTS DEFINED
#define LAMA_ASSERT(exp, msg)                                                           \
    if (!(exp))                                                                         \
    {                                                                                   \
        std::ostringstream errorStr;                                                    \
        errorStr << "Assertion Failed in line " << __LINE__;                            \
        errorStr << " of file " << __FILE__ << "\n";                                    \
        errorStr << "    Message: " << msg << "\n";                                     \
        lama::Exception::addCallStack( errorStr );                                      \
        throw lama::Exception( errorStr.str() );                                        \
    }

#define LAMA_ASSERT_EQUAL(exp1, exp2)                                                   \
    LAMA_ASSERT(exp1 == exp2, #exp1 " = " << exp1                                   \
                << " must be equal to " #exp2 " = " << exp2 )

#endif // LAMA_CHECK_ASSERTS
#if defined(LAMA_ASSERT_LEVEL_OFF)

#define LAMA_ASSERT_ERROR(exp, msg)                                                  \
    {                                                                               \
        unused( exp );                                                              \
        if ( false )                                                                \
        {                                                                           \
            std::cout<<msg;                                                         \
        }                                                                           \
    }
#define LAMA_ASSERT_EQUAL_ERROR(exp1, exp2) unused( exp1, exp2 );
#define LAMA_ASSERT_DEBUG(exp, msg)                                                  \
    {                                                                               \
        unused( exp );                                                              \
        if ( false )                                                                \
        {                                                                           \
            std::cout<<msg;                                                         \
        }                                                                           \
    }
#define LAMA_ASSERT_EQUAL_DEBUG(exp1, exp2) unused(exp1, exp2 );

#elif defined(LAMA_ASSERT_LEVEL_ERROR)

/**
 * @brief LAMA_ASSERT_ERROR(exp, msg) used for assertions that should cause an error
 */
#define LAMA_ASSERT_ERROR(exp, msg)  LAMA_ASSERT(exp, msg)
#define LAMA_ASSERT_DEBUG(exp, msg)  unused( exp );
#define LAMA_ASSERT_EQUAL_ERROR(exp1, exp2)  LAMA_ASSERT_EQUAL(exp1, exp2)
#define LAMA_ASSERT_EQUAL_DEBUG(exp1, exp2)  unused( exp1, exp2 );

#elif defined(LAMA_ASSERT_LEVEL_DEBUG)

/**
 * @brief LAMA_ASSERT_DEBUG(exp, msg) used in a debug mode for additional tests.
 */
#define LAMA_ASSERT_DEBUG(exp, msg)  LAMA_ASSERT(exp, msg)
#define LAMA_ASSERT_ERROR(exp, msg)  LAMA_ASSERT(exp, msg)

#define LAMA_ASSERT_EQUAL_ERROR(exp1, exp2)  LAMA_ASSERT_EQUAL(exp1, exp2)
#define LAMA_ASSERT_EQUAL_DEBUG(exp1, exp2)  LAMA_ASSERT_EQUAL(exp1, exp2)

#else

#pragma message( "Please specify LAMA_ASSERT_LEVEL_xxx with xxx = DEBUG, ERROR, or OFF" )
#pragma message( "Will use LAMA_ASSERT_LEVEL_ERROR by default." )

/**
 * @brief LAMA_ASSERT_ERROR(exp, msg) used for assertions that should cause an error
 */
#define LAMA_ASSERT_ERROR(exp, msg)  LAMA_ASSERT(exp, msg)
#define LAMA_ASSERT_DEBUG(exp, msg)  unused( exp );
#define LAMA_ASSERT_EQUAL_ERROR(exp1, exp2)  LAMA_ASSERT_EQUAL(exp1, exp2)
#define LAMA_ASSERT_EQUAL_DEBUG(exp1, exp2)  unused( exp1, exp2 );

#endif

#endif // LAMA_ASSERT_HPP_
