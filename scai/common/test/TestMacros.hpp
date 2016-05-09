/**
 * @file scai/common/test/TestMacros.hpp
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
 * @brief Additional Macros used for testing of LAMA with Boost Test.
 * @author Jiri Kraus
 * @date 06.04.2011
 * @since 1.0.0
 */

#pragma once

// local library
#include <scai/common/Printable.hpp>
#include <scai/common/exception/UnsupportedException.hpp>
#include <scai/common/SCAITypes.hpp>
#include <scai/common/ScalarType.hpp>
#include <scai/common/TypeTraits.hpp>
#include <scai/common/Math.hpp>

#include <scai/common/test/Configuration.hpp>

#include <scai/common/Math.hpp>
#include <scai/common/TypeTraits.hpp>

#include <scai/common/boost_check/version.hpp>

// boost
#include <boost/assign/list_of.hpp>
#include <boost/mpl/list.hpp>

/* -------------------------------------------------------------------------------- */
/*  arithmetic test types                                                           */
/* -------------------------------------------------------------------------------- */

typedef boost::mpl::list<SCAI_ARITHMETIC_HOST> scai_arithmetic_test_types;

/* -------------------------------------------------------------------------------- */
/*  array test types                                                                */
/* -------------------------------------------------------------------------------- */

typedef boost::mpl::list<SCAI_ARITHMETIC_ARRAY_HOST> scai_array_test_types;

/* -------------------------------------------------------------------------------- */
/*  Define one type used for test where only one ARITHMETIC TYPE is needed          */
/* -------------------------------------------------------------------------------- */

#define SCAI_TEST_TYPE SCAI_COMMON_FIRST_ENTRY( SCAI_ARITHMETIC_HOST )

/* -------------------------------------------------------------------------------- */

/*
 * @brief HelperMacro SCAI_CHECK_CLOSE( ValueType, x, y, tolerance )
 *
 * Extended macro BOOST_CHECK_CLOSE( x, y, tolerance) from Boost.Test.
 * as it does not work for complex types
 *
 * @param ValueType     value type to be used for check
 * @param x             Scalar
 * @param y             Scalar
 * @param percent_eps   Epsilon[%]
 *
 */

#define SCAI_CHECK_CLOSE( x, y, tolerance )                                                             \
    {                                                                                                   \
        BOOST_CHECK_CLOSE( scai::common::Math::real( x ), scai::common::Math::real( y ), tolerance );   \
        BOOST_CHECK_CLOSE( scai::common::Math::imag( x ), scai::common::Math::imag( y ), tolerance ) ;  \
    }

/* -------------------------------------------------------------------------------- */
/*                                                                                  */
/*      SCAI_CHECK_THROW   ->  replaces BOOST_CHECK_THROW                           */
/*                                                                                  */
/* -------------------------------------------------------------------------------- */

#if defined(SCAI_ASSERT_LEVEL_OFF)

// Do not check for exception if ASSERTions are switched off

#define SCAI_CHECK_THROW( stmt, exception )

#else

// use Boost macro for throwing exception

#define SCAI_CHECK_THROW( stmt, exception )    \
    BOOST_CHECK_THROW( stmt, exception )

#endif

/*
 * @brief TestMacro SCAI_COMMON_WRITEAT_TEST( printable )
 *
 * This macro checks if a output will be created by writing an object
 * into a stream. The length of this output must be greater than 0.
 * This object must be inherited from class Printable.
 *
 * @param printable     object of type printable
 */

#define SCAI_COMMON_WRITEAT_TEST( printable )              \
    {                                                      \
        std::stringstream mStream;                         \
        mStream << printable;                              \
        std::string mString = mStream.str();               \
        BOOST_CHECK( mString.length() > 0 );               \
    }

