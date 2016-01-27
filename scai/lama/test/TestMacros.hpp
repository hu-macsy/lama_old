/**
 * @file scai/lama/test/TestMacros.hpp
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

#include <scai/common/test/TestMacros.hpp>
#include <scai/hmemo/test/TestMacros.hpp>
#include <scai/kregistry/test/TestMacros.hpp>

#include <scai/lama/Scalar.hpp>

/*
 * Redefinition of SCAI_CHECK_CLOSE
 * works even if ValueType is unknown
 */

#ifdef SCAI_CHECK_CLOSE

#undef 	SCAI_CHECK_CLOSE

#endif

#ifdef SCAI_COMPLEX_SUPPORTED

#define SCAI_CHECK_CLOSE( x, y, tolerance )                                  \
    {                                                                        \
        scai::lama::Scalar xScalar = scai::lama::Scalar( x );                \
        scai::lama::Scalar yScalar = scai::lama::Scalar( y );                \
        ComplexDouble xVal = xScalar.getValue<ComplexDouble>();              \
        ComplexDouble yVal = yScalar.getValue<ComplexDouble>();              \
        BOOST_CHECK_CLOSE( xVal.real(), yVal.real(), tolerance );            \
        BOOST_CHECK_CLOSE( xVal.imag(), yVal.imag(), tolerance );            \
    }

#else

#define SCAI_CHECK_CLOSE( x, y, tolerance )                         \
    {                                                               \
        scai::lama::Scalar xScalar = scai::lama::Scalar( x );       \
        scai::lama::Scalar yScalar = scai::lama::Scalar( y );       \
        double xVal = xScalar.getValue<double>();                   \
        double yVal = yScalar.getValue<double>();                   \
        BOOST_CHECK_CLOSE( xVal.real(), yVal.real(), tolerance );   \
        BOOST_CHECK_CLOSE( xVal.imag(), yVal.imag(), tolerance );   \
    }

#endif

#define SCAI_CHECK_SMALL( x, ValueType, eps )                   \
        BOOST_CHECK_SMALL( x, static_cast<ValueType>(eps) );    \

#define SCAI_CHECK_SMALL_EPS( x, ValueType )                  \
    SCAI_CHECK_SMALL( x, ValueType, eps<ValueType> () )

/*
 * @brief HelperMacro SCAI_CHECK_SCALAR_SMALL( x, ValueType, eps )
 *
 * Extended Macro BOOST_CHECK_SMALL( x, eps ) from Boost.Test for
 * Scalar class of LAMA. Transforms Scalar x into ValueType,
 * and calls BOOST_CHECK_SMALL with arguments of type ValueType.
 *
 * @param x             Scalar
 * @param ValueType     type of Scalar x used for test
 * @param eps           Epsilon
 *
 * Static cast is used to convert eps to the right ValueType.
 */

/*#define SCAI_CHECK_SCALAR_SMALL( x, ValueType, eps )                     \
    {                                                                    \
        ValueType xHelper = (x).getValue<ValueType >();                  \
        BOOST_CHECK_SMALL( xHelper, static_cast<ValueType >( eps ) );    \
    }*/

#define SCAI_CHECK_SCALAR_SMALL( x, ValueType, eps )                     \
        SCAI_CHECK_SMALL( (x).getValue<ValueType>(), ValueType, eps )    \

/*
 * @brief HelperMacro SCAI_CHECK_SCALAR_SMALL_EPS( x, ValueType )
 *
 * Same as SCAI_CHECK_SCALAR_SMALL but with default eps value.
 *
 * @param x             Scalar
 * @param ValueType     type of Scalar to be used for test
 */

#define SCAI_CHECK_SCALAR_SMALL_EPS( x, ValueType )                  \
    SCAI_CHECK_SCALAR_SMALL( x, ValueType, eps<ValueType> () )

