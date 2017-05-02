/**
 * @file SparseKernelTest.cpp
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
 * @brief Contains tests for the SparseKernel interface to be tested on different devices
 * @author Thomas Brandes
 * @date 05.07.2013
 */

// boost
#include <boost/test/unit_test.hpp>

// others
#include <scai/hmemo.hpp>
#include <scai/kregistry.hpp>
#include <scai/utilskernel/LAMAKernel.hpp>
#include <scai/utilskernel/LArray.hpp>
#include <scai/sparsekernel/StencilKernelTrait.hpp>
#include <scai/common/Settings.hpp>
#include <scai/common/test/TestMacros.hpp>

#include <scai/hmemo/test/ContextFix.hpp>

/*--------------------------------------------------------------------- */

using namespace scai;
using namespace hmemo;
using namespace sparsekernel;
using namespace utilskernel;
using common::TypeTraits;

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( SparseKernelTest )

/* --------------------------------------------------------------------- */

SCAI_LOG_DEF_LOGGER( logger, "Test.SparseKernelTest" )

/* ------------------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( stencilLocalTest, ValueType, scai_numeric_test_types )
{
    ContextPtr testContext = ContextFix::testContext;
    ContextPtr hostContext = Context::getHostPtr();

    LAMAKernel<StencilKernelTrait::stencilLocalSizes> stencilLocalSizes;
    LAMAKernel<StencilKernelTrait::stencilLocalCSR<ValueType> > stencilLocalCSR;

    ContextPtr loc = testContext;
    stencilLocalSizes.getSupportedContext( loc, stencilLocalCSR );

    LArray<IndexType> csrIA( testContext );
    LArray<IndexType> csrJA( testContext );
    LArray<ValueType> csrValues( testContext );

    IndexType nDims = 1;
    IndexType gridSizes[] = { 100 };
    IndexType gridDistances[] = { 1 };
    IndexType nPoints = 3;
    int stencilNodes[] = { -1, 0, 1 };
    // ValueType stencilValues[] = { -1, 2, -1 };

    {
        IndexType n = gridSizes[0];   // number of grid points
        WriteOnlyAccess<IndexType> wIA( csrIA, loc, n + 1 );
        stencilLocalSizes[loc]( wIA.get(), nDims, gridSizes, gridDistances, nPoints, stencilNodes );
        wIA.resize( n );
    }

    BOOST_CHECK_EQUAL( 3, csrIA[1] );   // inner point: all 3 stencil points are valid
    BOOST_CHECK_EQUAL( 2, csrIA[0] );   // left boundary, only 2 points
    BOOST_CHECK_EQUAL( 2, csrIA[99] );  // right boundary, only 2 points

    // now build the CSR offsets

    IndexType nnz = HArrayUtils::scan1( csrIA );

    {
        WriteOnlyAccess<IndexType> wJA( csrJA, loc, nnz );
        WriteOnlyAccess<ValueType> wValues( csrValues, loc, nnz );
    }
}

/* ------------------------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END()
