/**
 * @file LAMAAssert.hpp
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
 * @brief Contains error checking macros like LAMA_ASSERT that can be enabled
 *        or disabled at runtime
 *
 * @author Jiri Kraus
 * @date 02.03.2011
 */

#pragma once

#include <lama/exception/Exception.hpp>
#include <common/Exception.hpp>

#include <sstream>
#include <cstdio>
#include <iostream>

// Little help for the Eclipse Parser

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

#define LAMA_ASSERT(exp, msg) unused( exp );

#define LAMA_ASSERT_EQUAL(exp1, exp2) unused( exp1, exp2 );

#else // LAMACHECKASSERTS DEFINED

#define LAMA_ASSERT(exp, msg)                                                           \
    COMMON_ASSERT( exp, msg )

#define LAMA_ASSERT_EQUAL(exp1, exp2)                                                   \
    COMMON_ASSERT_EQUAL( exp1, exp2, "LAMA equality error" ) 

#endif // LAMA_CHECK_ASSERTS

#if defined( LAMA_ASSERT_LEVEL_OFF )

#define LAMA_ASSERT_ERROR(exp, msg)                                                 \
    {                                                                               \
        unused( exp );                                                              \
        if ( false )                                                                \
        {                                                                           \
            std::cout<<msg;                                                         \
        }                                                                           \
    }

#define LAMA_ASSERT_EQUAL_ERROR( exp1, exp2 )                                       \
    unused( exp1, exp2 );

#define LAMA_ASSERT_DEBUG(exp, msg)                                                 \
    {                                                                               \
        unused( exp );                                                              \
        if ( false )                                                                \
        {                                                                           \
            std::cout<<msg;                                                         \
        }                                                                           \
    }

#define LAMA_ASSERT_EQUAL_DEBUG( exp1, exp2 )                                      \
    unused( exp1, exp2 );

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
