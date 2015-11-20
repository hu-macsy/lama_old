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

/*
 * Redefinition of SCAI_CHECK_CLOSE
 * works even if ValueType is unknown
 */

#ifdef SCAI_CHECK_CLOSE

#undef 	SCAI_CHECK_CLOSE

#define SCAI_CHECK_CLOSE( x, y, tolerance )                         \
    {                                                               \
        Scalar xScalar = Scalar( x );                               \
        Scalar yScalar = Scalar( y );                               \
        ComplexDouble xVal = xScalar.getValue<ComplexDouble>();     \
        ComplexDouble yVal = yScalar.getValue<ComplexDouble>();     \
        BOOST_CHECK_CLOSE( xVal.real(), yVal.real(), tolerance );   \
        BOOST_CHECK_CLOSE( xVal.imag(), yVal.imag(), tolerance );   \
    }

#endif
