/**
 * @file scai/common/test/TestMacros.hpp
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
 * @brief Additional Macros used for testing of LAMA with Boost Test.
 * @author Jiri Kraus
 * @date 06.04.2011
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

#include <scai/common/boost_check.hpp>

// boost
#include <boost/assign/list_of.hpp>
#include <boost/mpl/list.hpp>

/* -------------------------------------------------------------------------------- */
/*  arithmetic test types                                                           */
/* -------------------------------------------------------------------------------- */

typedef boost::mpl::list<SCAI_NUMERIC_TYPES_HOST> scai_numeric_test_types;

/* -------------------------------------------------------------------------------- */
/*  array test types                                                                */
/* -------------------------------------------------------------------------------- */

typedef boost::mpl::list<SCAI_ARRAY_TYPES_HOST> scai_array_test_types;

/* -------------------------------------------------------------------------------- */
/*  fft test types                                                                */
/* -------------------------------------------------------------------------------- */

typedef boost::mpl::list<SCAI_FFT_TYPES_HOST> scai_fft_test_types;

/* -------------------------------------------------------------------------------- */
/*  Define one type used for test where only one ARITHMETIC TYPE is needed          */
/* -------------------------------------------------------------------------------- */

#define SCAI_TEST_TYPE SCAI_COMMON_FIRST_ARG( SCAI_NUMERIC_TYPES_HOST )

/* -------------------------------------------------------------------------------- */

/*
 * @brief HelperMacro SCAI_CHECK_CLOSE( x, y, tolerance )
 *
 * Extends macro BOOST_CHECK_CLOSE( x, y, tolerance) from Boost.Test.
 * as it does not work for complex types
 *
 * @param x             Scalar
 * @param y             Scalar
 * @param tolerance     Epsilon[%]
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

