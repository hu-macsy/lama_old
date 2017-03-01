/**
 * @file SparseVectorTest.cpp
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
 * @brief Contains specific tests for class SparseVector (mainly constructors)
 * @author Thomas Brandes
 * @date 17.01.2017
 */

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <scai/lama/test/TestMacros.hpp>
#include <scai/lama/test/matrix/Matrices.hpp>

#include <scai/dmemo/test/TestDistributions.hpp>
#include <scai/dmemo/BlockDistribution.hpp>
#include <scai/dmemo/CyclicDistribution.hpp>

#include <scai/lama/SparseVector.hpp>
#include <scai/lama/matrix/SparseMatrix.hpp>
#include <scai/lama/matrix/CSRSparseMatrix.hpp>
#include <scai/lama/matutils/MatrixCreator.hpp>

using namespace scai;
using namespace lama;

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( SparseVectorTest )

/* --------------------------------------------------------------------- */

SCAI_LOG_DEF_LOGGER( logger, "Test.SparseVectorTest" )

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( cTorTest, ValueType, scai_numeric_test_types )
{
    // replicated sparse vector of a certain size only with zero elements

    IndexType n = 4;

    SparseVector<ValueType> v( n );

    ValueType zero = 0;

    BOOST_CHECK_EQUAL( n, v.size() );

    for ( IndexType i = 0; i < n; ++i )
    {
        BOOST_CHECK_EQUAL( v.getValue( i ), zero );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( consistencyTest )
{
    typedef RealType ValueType;

    IndexType n = 10;

    dmemo::CommunicatorPtr comm = dmemo::Communicator::getCommunicatorPtr();
    dmemo::DistributionPtr dist( new dmemo::BlockDistribution( n, comm ) );

    // create distributed sparse vector, initialized with 0

    SparseVector<ValueType> v( dist );

    BOOST_CHECK( v.isConsistent() );

    // usually it is not simple to make a vector inconsistent

    hmemo::HArray<IndexType>& indexes = const_cast<hmemo::HArray<IndexType>&>( v.getNonZeroIndexes() );

    // make it inconsistent on one processor only, isConsistent can deal with it

    if ( comm->getRank() == ( comm->getSize() / 2 ) )
    {
        indexes.resize( 1 );
    }

    // consistency check fails on all processors

    BOOST_CHECK( ! v.isConsistent() );

    hmemo::HArray<ValueType>& values = const_cast<hmemo::HArray<ValueType>&>( v.getNonZeroValues() );

    if ( dist->getLocalSize() > 0 )
    {
        indexes.init( IndexType( 0 ), 1 );
        values.init( ValueType( 3 ), 1 );
    }
    else
    {
        indexes.clear();
        values.clear();
    }

    BOOST_CHECK( v.isConsistent() );

    if ( dist->getLocalSize() > 0 )
    {
        indexes.init( 1, IndexType( dist->getLocalSize() ) );
    }

    BOOST_CHECK( ! v.isConsistent() );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( CopyConstructorTest )
{
    // Note: it is sufficient to consider one value type

    typedef RealType ValueType;

    IndexType n = 20;   // let it really small, this test is very inefficient

    utilskernel::LArray<ValueType> data;

    std::srand( 13151 );   // This test only works if all processors have same random numbers

    const float fillRate = 0.2;   // makes it worth to be sparse

    utilskernel::HArrayUtils::setRandom( data, n, fillRate );

    dmemo::TestDistributions dists( n );

    for ( size_t i = 0; i < dists.size(); ++i )
    {
        utilskernel::LArray<ValueType> data1( n );

        dmemo::DistributionPtr dist = dists[i];

        DenseVector<ValueType> denseV( data );

        denseV.redistribute( dist );

        SparseVector<ValueType> sparseV( denseV );

        // get each value from the distributed sparse vector

        for ( IndexType k = 0; k < n; ++k )
        {
            data1[k] = sparseV.getValue( k ).getValue<ValueType>();
        }

        BOOST_CHECK_EQUAL( 0, data.maxDiffNorm( data1 ) );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( SparseConstructorTest )
{
    typedef RealType ValueType;

    ValueType values_raw[] = { 3, 4, 3, 1, 2 };
    IndexType indexes_raw[] = { 7, 3, 11, 13, 5 };

    IndexType nnz = sizeof( values_raw ) / sizeof( ValueType );
    IndexType n   = 20;

    hmemo::HArray<ValueType> values( nnz, values_raw );
    hmemo::HArray<IndexType> indexes( nnz, indexes_raw );

    dmemo::DistributionPtr dist( new dmemo::NoDistribution( n ) );

    // The constructor copies the arrays and sorts the entries

    SparseVector<ValueType> s( indexes, values, dist );

    const hmemo::HArray<IndexType>& spIndexes = s.getNonZeroIndexes();

    BOOST_CHECK( utilskernel::HArrayUtils::isSorted( spIndexes, utilskernel::binary::LT ) );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( RedistributeTest )
{
    // Note: it is sufficient to consider one value type

    typedef RealType ValueType;

    IndexType n = 20;   // let it really small, this test is very inefficient

    utilskernel::LArray<ValueType> data;

    std::srand( 13151 );   // This test only works if all processors have same random numbers

    const float fillRate = 0.2;   // makes it worth to be sparse

    utilskernel::HArrayUtils::setRandom( data, n, fillRate );

    dmemo::TestDistributions dists( n );

    for ( size_t i = 0; i < dists.size(); ++i )
    {
        utilskernel::LArray<ValueType> data1( n );

        dmemo::DistributionPtr dist = dists[i];

        SparseVector<ValueType> sparseV( data );

        sparseV.redistribute( dist );

        // get each value from the distributed sparse vector

        for ( IndexType k = 0; k < n; ++k )
        {
            data1[k] = sparseV.getValue( k ).getValue<ValueType>();
        }

        BOOST_CHECK_EQUAL( 0, data.maxDiffNorm( data1 ) );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();
