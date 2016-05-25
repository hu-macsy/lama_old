/**
 * @file P_NormTest.cpp
 *
 * @license
 * Copyright (c) 2009-2016
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the Library of Accelerated Math Applications (LAMA).
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
 * @endlicense
 *
 * @brief Test norm for distributed vectors.
 * @author Thomas Brandes
 * @date 10.05.2013
 */

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <scai/lama/DenseVector.hpp>
#include <scai/dmemo/BlockDistribution.hpp>

#include <scai/lama/test/EquationHelper.hpp>
#include <scai/lama/test/TestMacros.hpp>

using namespace scai::lama;
using namespace scai::hmemo;
using namespace scai::dmemo;

typedef boost::mpl::list<float, double> test_types;

/* --------------------------------------------------------------------- */

BOOST_TEST_SUITE( NormTest )

SCAI_LOG_DEF_LOGGER( logger, "Test.P_NormTest" );

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( Norm, ValueType, test_types )
{
    ContextPtr context = Context::getContextPtr();
    CommunicatorPtr comm = Communicator::getCommunicator();

        for ( IndexType size = 0; size < 4; ++size )
        {
            DistributionPtr dist( new BlockDistribution( size, comm ) );
            const ValueType VAL = 1.0;
            DenseVector<ValueType> repVector( size, VAL );
            DenseVector<ValueType> distVector;
            distVector.setContextPtr( context );
            distVector = repVector;
            distVector.redistribute( dist );
            Scalar l1norm = distVector.l1Norm();
            ValueType expectedL1Norm = static_cast<ValueType>( size );
            BOOST_CHECK_CLOSE( l1norm.getValue<ValueType>(), expectedL1Norm, 1 );
            Scalar l2norm = distVector.l2Norm();
            ValueType expectedL2Norm = static_cast<ValueType>( size );
            expectedL2Norm = std::sqrt( expectedL2Norm );
            BOOST_CHECK_CLOSE( l2norm.getValue<ValueType>(), expectedL2Norm, 1 );
            Scalar maxNorm = distVector.maxNorm();
            ValueType expectedMaxNorm = static_cast<ValueType>( VAL );

            if ( size == 0 )
            {
                expectedMaxNorm = static_cast<ValueType>( 0 );
            }

            BOOST_CHECK_CLOSE( maxNorm.getValue<ValueType>(), expectedMaxNorm, 1 );
        }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();
