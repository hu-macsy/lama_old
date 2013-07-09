/**
 * @file BLAS1Test.cpp
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
 * @brief Contains tests for the class CUDABLAS1 and OpenMPBLAS1
 * @author: Bea Hornef
 * @date 5.7.2013
 * @since 1.0.0
 **/

// math for sqrt
#include <cmath>

// boost
#include <boost/test/unit_test.hpp>

// others
#include <lama/ContextAccess.hpp>
#include <lama/HostReadAccess.hpp>
#include <lama/LAMAArray.hpp>
#include <lama/LAMAInterface.hpp>
#include <lama/ReadAccess.hpp>
#include <lama/Scalar.hpp>
#include <lama/WriteAccess.hpp>

#include <test/TestMacros.hpp>

using namespace boost;
using namespace lama;

namespace lama
{
namespace BLAS1Test
{

template<typename ValueType>
void nrm2Test( ContextPtr loc )
{
    LAMA_INTERFACE_FN_T( nrm2, loc, BLAS, BLAS1, ValueType );

    {
        ValueType values[] =
        { 1, 2, 3, 4, 5, 6 };
        const IndexType nValues = sizeof( values ) / sizeof( ValueType );
        const IndexType incX1 = 1;
        const IndexType incX2 = 2;
        const ValueType result1 = 91.0;
        const ValueType result2 = 35.0;

        LAMAArray<ValueType> AValues( nValues, values );

        {
            LAMA_CONTEXT_ACCESS( loc );

            // std::cout << "test 1 (incX = 1)" << std::endl;
            ReadAccess<ValueType> rAValues( AValues, loc );
        	ValueType euclideanNorm = nrm2( nValues / incX1, rAValues.get(), incX1, NULL );
            BOOST_CHECK_CLOSE( euclideanNorm, ::sqrt(result1), 1e-4 );

            // std::cout << "test 2 (incX = 2)" << std::endl;
        	euclideanNorm = nrm2( nValues / incX2, rAValues.get(), incX2, NULL );
        	BOOST_CHECK_CLOSE( euclideanNorm, ::sqrt(result2), 1e-4 );
        }
    }
}

} // namespace BLAS1Test
} // namespace lama

/* ------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_SUITE( BLAS1Test );

LAMA_LOG_DEF_LOGGER( logger, "Test.BLAS1Test" );

LAMA_AUTO_TEST_CASE_T( nrm2Test, BLAS1Test );

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_SUITE_END();
