/**
 * @file VectorTest.cpp
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
 * @brief Contains generic tests for Vector objects.
 * @author Thomas Brandes
 * @date 27.07.2016
 */

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <scai/lama/test/TestMacros.hpp>
#include <scai/lama/test/TestVectors.hpp>

#include <scai/dmemo/test/TestDistributions.hpp>
#include <scai/dmemo/BlockDistribution.hpp>

#include <scai/lama/DenseVector.hpp>
#include <scai/lama/matrix/CSRSparseMatrix.hpp>

using namespace scai;
using namespace lama;

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( VectorTest )

/* --------------------------------------------------------------------- */

SCAI_LOG_DEF_LOGGER( logger, "Test.VectorTest" )

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( WriteTest )
{
    const IndexType n = 13;

    TestVectors vectors;

    for ( size_t i = 0; i < vectors.size(); ++i )
    {
        VectorPtr v = vectors[i];
        v->allocate( n );

        std::ostringstream out1;
        out1 << *v;   // same as v1->writeAt( out1 );
        BOOST_CHECK( out1.str().length() > 0 );

        std::ostringstream out2;
        v->Vector::writeAt( out2 );
        BOOST_CHECK( out1.str().length() > 0 );

        BOOST_CHECK( out2.str() != out1.str() );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( AllocateTest )
{
    dmemo::CommunicatorPtr comm( dmemo::Communicator::getCommunicatorPtr() );

    const IndexType n = 13;

    TestVectors vectors;

    dmemo::DistributionPtr dist( new dmemo::BlockDistribution( n, comm ) );
    dmemo::DistributionPtr dist1( new dmemo::BlockDistribution( n + 1, comm ) );
    dmemo::DistributionPtr repDist1( new dmemo::NoDistribution( n + 1 ) );

    for ( size_t i = 0; i < vectors.size(); ++i )
    {
        VectorPtr v = vectors[i];

        size_t size0 = v->getMemoryUsage();

        v->allocate( n );

        size_t size1 = v->getMemoryUsage();

        BOOST_CHECK( size1 >= size0 + n * common::typeSize( v->getValueType() ) );

        *v = 1;
        v->redistribute( dist );
        BOOST_CHECK_EQUAL( ( *v )( n - 1 ), Scalar( 1 ) );

        BOOST_CHECK_THROW( 
        {
           v->redistribute( dist1 );
        }, common::Exception );

        v->allocate( dist1 );
        *v = 2;
        v->redistribute( repDist1 );

        BOOST_CHECK_EQUAL( ( *v )( n - 1 ), Scalar( 2 ) );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( InvertTest )
{   
    dmemo::CommunicatorPtr comm( dmemo::Communicator::getCommunicatorPtr() );
    
    const IndexType n = 13;
    
    TestVectors vectors;
    
    dmemo::DistributionPtr dist( new dmemo::BlockDistribution( n, comm ) ); 
    
    for ( size_t i = 0; i < vectors.size(); ++i )
    {   
        VectorPtr v = vectors[i];
        
        if ( ! common::isNumeric( v->getValueType() ) )
        {
            continue;   // this test does not work for int, uint, ....
        }

        v->allocate( dist );

        *v = 4;
 
        v->invert();
 
        Scalar s = ( *v )( n / 2 );

        // s should be 2, but might not be exact

        BOOST_CHECK( ( s - Scalar( 0.25) ) < 0.00001 );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( assign_S_VV_Test )
{
    const IndexType n = 13;

    TestVectors vectors;

    dmemo::TestDistributions dists(n);

    for ( size_t i = 0; i < vectors.size(); ++i )
    {
        VectorPtr v1 = vectors[i];

        if ( ! common::isNumeric( v1->getValueType() ) )
        {
            continue;   // this test does not work for int, uint, ....
        }

        for ( size_t j = 0; j < dists.size(); ++j )
        {
            dmemo::DistributionPtr dist = dists[j];

            v1->allocate( dist );
            *v1 = 3;
            VectorPtr v2( v1->copy() );
            *v2 = 5;
            VectorPtr v3( v1->newVector() );
            VectorPtr v4( v1->newVector() );
            *v3 = *v1 * *v2;
            *v4 = 3 * *v1 * *v2;
            *v4 -= 2 * *v3;
 
            // Now v3 and v4 must be equal

            *v3 -= *v4;

            BOOST_CHECK( v3->maxNorm() < Scalar( 1e-4 ) );
        }
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( assign_MV_Test )
{
    const IndexType n = 13;

    TestVectors vectors;

    dmemo::TestDistributions dists(n);

    for ( size_t i = 0; i < vectors.size(); ++i )
    {
        VectorPtr v1 = vectors[i];

        if ( ! common::isNumeric( v1->getValueType() ) )
        {
            continue;   // this test does not work for int, uint, ....
        }

        for ( size_t j = 0; j < dists.size(); ++j )
        {
            dmemo::DistributionPtr dist = dists[j];

            v1->setSequence( 3, 1, dist );

            MatrixPtr m( Matrix::getMatrix( Matrix::CSR, v1->getValueType() ));
            m->setIdentity( dist );

            VectorPtr v2( v1->newVector() );
 
            *v2 = *m * *v1;

            // Now v1 and v2 must be equal

            *v2 -= *v1;

            BOOST_CHECK( v2->maxNorm() < Scalar( 1e-4 ) );

            *v2 = 2 * *m * *v1;

            // Now v1 and v2 must be equal

            *v2 -= 2 * *v1;

            BOOST_CHECK( v2->maxNorm() < Scalar( 1e-4 ) );

            *v2 = *m * *v1 - *v1;
            
            BOOST_CHECK( v2->maxNorm() < Scalar( 1e-4 ) );

            *v2 = 2 * *m * *v1 - 2 * *v1;
            
            BOOST_CHECK( v2->maxNorm() < Scalar( 1e-4 ) );
        }
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( assign_VM_Test )
{
    return;  

    // This test fails sometimes with 5 or 6 processors
    // valgrind shows memory problems during MPI gather
    // TODO: Lauretta  

    const IndexType n = 13;

    TestVectors vectors;

    dmemo::TestDistributions dists(n);

    for ( size_t i = 0; i < vectors.size(); ++i )
    {
        VectorPtr v1 = vectors[i];

        if ( ! common::isNumeric( v1->getValueType() ) )
        {
            continue;   // this test does not work for int, uint, ....
        }

        for ( size_t j = 0; j < dists.size(); ++j )
        {   
            dmemo::DistributionPtr dist = dists[j];
            
            v1->allocate( dist );
            *v1 = 3;
            
            MatrixPtr m( Matrix::getMatrix( Matrix::CSR, v1->getValueType() ));
            m->setIdentity( dist );
            m->setCommunicationKind( Matrix::ASYNCHRONOUS );
            
            // SCAI_LOG_ERROR( logger, "vectorMultMatrix, v = " << *v1 << ", m = " << *m )

            VectorPtr v2( v1->newVector() );

            *v2 = ( *v1 ) * ( *m );

            // Now v1 and v2 must be equal

            *v2 -= *v1;
            
            Scalar s = v2->maxNorm();

            BOOST_CHECK( s < Scalar( 1e-4 ) );

            *v2 = 2 * *v1 * *m;
            
            // Now 2 * v1 and v2 must be equal
            
            *v2 -= 2 * *v1;
            
            BOOST_CHECK( v2->maxNorm() < Scalar( 1e-4 ) );
            
            *v2 = *v1 * *m - *v1;
            
            BOOST_CHECK( v2->maxNorm() < Scalar( 1e-4 ) );
            
            *v2 = 2 * *v1 * *m - 2 * *v1;
            
            BOOST_CHECK( v2->maxNorm() < Scalar( 1e-4 ) );
        }
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( GetDenseVectorTest, ValueType, scai_numeric_test_types )
{
    hmemo::ContextPtr ctx = hmemo::Context::getContextPtr();
    dmemo::CommunicatorPtr comm = dmemo::Communicator::getCommunicatorPtr();

    IndexType n = 111;

    dmemo::DistributionPtr dist( new dmemo::BlockDistribution( n, comm ) );

    common::scalar::ScalarType stype = common::TypeTraits<ValueType>::stype;

    VectorPtr v( Vector::getDenseVector( stype, dist, ctx ) );

    std::string format = Vector::kind2Str( v->getVectorKind() );

    BOOST_CHECK_EQUAL( "DENSE", format );
    BOOST_CHECK_EQUAL( Vector::str2Kind( "DENSE" ), v->getVectorKind() );
    
    DenseVector<ValueType>* denseV = dynamic_cast<DenseVector<ValueType>* >( v.get() );

    BOOST_REQUIRE( denseV != NULL );

    BOOST_CHECK_EQUAL( denseV->getContextPtr(), ctx );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( dotProductTest )
{   
    // This test does not verify the correctness of the dot product
    // but that it is the same, either computed replicated or distributed.

    using namespace hmemo;

    const IndexType n = 13;
    
    TestVectors vectors;
    
    dmemo::TestDistributions dists(n);
    
    std::srand( 1311 );   // same random numbers on all processors

    for ( size_t i = 0; i < vectors.size(); ++i )
    {   
        VectorPtr v1 = vectors[i];

        if ( ! common::isNumeric( v1->getValueType() ) )
        {
            continue;   // this test does not work for int, uint, ....
        }

        VectorPtr v2( v1->newVector() );

        // generate two arrays of same value type with random numbers

        common::shared_ptr<_HArray> data1( _HArray::create( v1->getValueType() ) );
        common::shared_ptr<_HArray> data2( _HArray::create( v1->getValueType() ) );

        utilskernel::HArrayUtils::setRandom( *data1, n );
        utilskernel::HArrayUtils::setRandom( *data2, n );

        v1->assign( *data1 );
        v2->assign( *data2 );

        Scalar dotp = v1->dotProduct( *v2 );  // replicated computation
 
        for ( size_t j = 0; j < dists.size(); ++j )
        {   
            // now compute the dot product with distributed vectors

            dmemo::DistributionPtr dist = dists[j];
            
            v1->redistribute( dist );
            v2->redistribute( dist );
            
            Scalar distDotp = v1->dotProduct( *v2 );

            // we cannot check for equality due to different rounding errors

            SCAI_CHECK_CLOSE( dotp, distDotp, 0.001 );
        }
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( scaleTest )
{   
    using namespace hmemo;

    const IndexType n = 13;
    
    TestVectors vectors;
    
    dmemo::TestDistributions dists(n);
    
    std::srand( 1311 );   // same random numbers on all processors

    for ( size_t i = 0; i < vectors.size(); ++i )
    {   
        VectorPtr v1 = vectors[i];

        if ( ! common::isNumeric( v1->getValueType() ) )
        {
            continue;   // this test does not work for int, uint, ....
        }

        common::shared_ptr<_HArray> data1( _HArray::create( v1->getValueType() ) );
        utilskernel::HArrayUtils::setRandom( *data1, n );

        v1->assign( *data1 );

        for ( size_t j = 0; j < dists.size(); ++j )
        {
            v1->redistribute( dists[j] );

            // CSRSparseMatrix<ValueType> m;  

            MatrixPtr m( Matrix::getMatrix( Matrix::CSR, v1->getValueType() ) );

            m->setIdentity( dists[j] );
            m->setDiagonal( *v1 );

            VectorPtr v2( v1->newVector() );

            *v2 = *v1;
            *v2 += *v1;
            *v2 -= *v1;
            *v2 *= 2;
            *v2 /= 2;

            *v2 *= *v2;   // is v1 * v1

            VectorPtr v3( v1->newVector() );

            *v3 = 2 * ( *m ) * ( *v1 ) - ( *v2 );   // is v1 * v1
 
            // ToDO: Lauretta: v3 -= v1 * m  

            *v3 -= *v2;

            BOOST_CHECK( v3->maxNorm() < Scalar( 0.0001 ) );
        }
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();
