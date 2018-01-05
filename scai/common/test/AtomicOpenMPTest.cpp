/**
 * @file atomicOpenMPTest.cpp
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
 * @brief Test routines for atomic openmp wrapper
 * @author Lauretta Schubert
 * @date 30.03.2016
 */

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <scai/common/test/TestMacros.hpp>
#include <scai/common/SCAITypes.hpp>

#include <scai/common/OpenMP.hpp>

using namespace scai;
using namespace common;

BOOST_AUTO_TEST_SUITE( AtomicOpenMPTest )

/* -------------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( atomicAddTest, ValueType, scai_numeric_test_types )
{
    int size = 100;
    ValueType globalResult = 0;
    #pragma omp parallel for

    for ( int i = 0; i < size; ++i )
    {
        ValueType localResult = i + static_cast<ValueType> ( 1 );
        atomicAdd( globalResult, localResult );
    }

    ValueType res = ( size * ( size + 1 ) ) / static_cast<ValueType> ( 2 );
    BOOST_CHECK_EQUAL( globalResult, res );
}

BOOST_AUTO_TEST_SUITE_END()
