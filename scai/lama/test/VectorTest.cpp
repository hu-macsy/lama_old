/**
 * @file VectorTest.cpp
 *
 * @license
 * Copyright (c) 2009-2016
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
 * @brief Contains generic tests for Vector (mainly constructors)
 * @author Thomas Brandes
 * @date 27.07.2016
 */

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <scai/lama/test/TestMacros.hpp>
#include <scai/lama/test/TestVectors.hpp>

#include <scai/dmemo/test/TestDistributions.hpp>

#include <scai/lama/DenseVector.hpp>
#include <scai/lama/matrix/CSRSparseMatrix.hpp>

using namespace scai;
using namespace lama;

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( VectorTest )

/* --------------------------------------------------------------------- */

SCAI_LOG_DEF_LOGGER( logger, "Test.VectorTest" )

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( vecExpAssignTest )
{
    const IndexType n = 13;

    TestVectors vectors;

    dmemo::TestDistributions dists(n);

    for ( size_t i = 0; i < vectors.size(); ++i )
    {
        VectorPtr v1 = vectors[i];

        for ( size_t j = 0; j < dists.size(); ++j )
        {
            dmemo::DistributionPtr dist = dists[i];

            v1->allocate( dist );
            *v1 = 3;

            SCAI_LOG_ERROR( logger, "run vec exp test for v = " << *v1 )

            VectorPtr v2( v1->copy() );
            *v2 = 5;

            VectorPtr v3( v1->newVector() );
            *v3 = *v1 + *v2;
            VectorPtr v4( v1->copy() );
            *v4 += *v2;
            *v3 -= 2 * *v4;
            BOOST_CHECK( v3->maxNorm() < Scalar( 1e-4 ) );
        }
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();
