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
#include <scai/lama/SparseVector.hpp>
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

        if ( v->getVectorKind() == Vector::DENSE )
        {
            // a dense vector allocates really memory

            BOOST_CHECK( size1 >= size0 + n * common::typeSize( v->getValueType() ) );
        }
        else if ( v->getVectorKind() == Vector::SPARSE )
        {
            // a sparse vector does not allocate here memory

            BOOST_CHECK_EQUAL( size1, size0 );
        }

        *v = 1;
        v->redistribute( dist );
        BOOST_CHECK_EQUAL( v->getValue( n - 1 ), Scalar( 1 ) );

        BOOST_CHECK_THROW(
        {
            v->redistribute( dist1 );
        }, common::Exception );

        v->allocate( dist1 );
        *v = 2;
 

        v->redistribute( repDist1 );

        continue;

        BOOST_CHECK_EQUAL( v->getValue( n - 1 ), Scalar( 2 ) );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( SetGetTest )
{
    dmemo::CommunicatorPtr comm( dmemo::Communicator::getCommunicatorPtr() );

    const IndexType n = 13;

    TestVectors vectors;

    for ( size_t i = 0; i < vectors.size(); ++i )
    {
        Vector& v = *vectors[i];

        v.allocate( n );

        v = 1;

        Scalar s = v[0];
        BOOST_CHECK_EQUAL( s, Scalar( 1 ) );

        v[n-2] = 9;
        v[1] = 7;

        BOOST_CHECK_THROW(
        {
            v[n] = 1;
        }, common::Exception );

        s = v[2];
        BOOST_CHECK_EQUAL( s, Scalar( 1 ) );
        s = v[1];
        BOOST_CHECK_EQUAL( s, Scalar( 7 ) );
        v[1] = 5;
        s = v[1];
        BOOST_CHECK_EQUAL( s, Scalar( 5 ) );
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

        BOOST_CHECK( ( s - Scalar( 0.25 ) ) < 0.00001 );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( ConjTest )
{
    dmemo::CommunicatorPtr comm( dmemo::Communicator::getCommunicatorPtr() );

    const IndexType n = 100;

    TestVectors vectors;

    dmemo::DistributionPtr dist( new dmemo::BlockDistribution( n, comm ) );

    for ( size_t i = 0; i < vectors.size(); ++i )
    {
        VectorPtr v = vectors[i];

        if ( ! common::isComplex( v->getValueType() ) )
        {
            continue;   // test only useful for complex numbers
        }

        std::srand( i + 3351 );  // same values on each processor

        float fillRate = 0.1f;

        v->setSparseRandom( dist, 0, fillRate, 1 );
        
        VectorPtr v1( v->copy() );

        v->conj();

        *v *= *v1;  // ( a + b i ) ( a - b i )
        
        Scalar s1 = v->sum();
        Scalar s2 = v1->dotProduct( *v1 );

        SCAI_LOG_DEBUG( logger, "sum( v * conj(v ) = " << s1 << ", dotProduct( v, v ) = " << s2 )

        BOOST_CHECK( abs( s1 - s2 ) < Scalar( 0.0001 ) );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( ExpLogTest )
{
    dmemo::CommunicatorPtr comm( dmemo::Communicator::getCommunicatorPtr() );

    const IndexType n = 100;

    TestVectors vectors;

    dmemo::DistributionPtr dist( new dmemo::BlockDistribution( n, comm ) );

    for ( size_t i = 0; i < vectors.size(); ++i )
    {
        Vector& v1 = *vectors[i];

        if ( ! common::isNumeric( v1.getValueType() ) )
        {
            continue;   // this test will fail for IndexType
        }

        float fillRate = 0.1f;

        Scalar zero = 0;

        v1.setSparseRandom( dist, zero, fillRate, 1 );

        v1 += 2;

        VectorPtr v2Ptr( v1.copy() );
        const Vector& v2 = *v2Ptr;

        v1.exp();
        v1.log();

        v1 -= v2;

        Scalar diff = v1.maxNorm();

        SCAI_LOG_DEBUG( logger, "v = " << v1 << ": maxNorm( log( exp ( v ) ) - v ) = " << diff )

        BOOST_CHECK( diff < Scalar( 0.0001 ) );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( SinCosTest )
{
    dmemo::CommunicatorPtr comm( dmemo::Communicator::getCommunicatorPtr() );

    const IndexType n = 100;

    TestVectors vectors;

    dmemo::DistributionPtr vectorDist( new dmemo::BlockDistribution( n, comm ) );

    for ( size_t i = 0; i < vectors.size(); ++i )
    {
        Vector& v1 = *vectors[i];

        if ( ! common::isNumeric( v1.getValueType() ) )
        {
            continue;   // this test would fail for IndexType
        }

        float fillRate = 0.1f;

        v1.setSparseRandom( vectorDist, 0, fillRate, 1 );

        VectorPtr v2Ptr( v1.copy() );
        Vector& v2 = *v2Ptr;

        // build:  sin(v1) * sin(v1) + cos(v2) * cos(v2) - 1, must all be 0

        v1.sin();
        v2.cos();

        // v1 = v1 * v1 - v2 * v2 - 1

        v1 *= v1;
        v2 *= v2;
        v1 += v2;
        v1 -= Scalar( 1 );

        Scalar diff = v1.maxNorm();

        BOOST_CHECK( diff < Scalar( 0.0001 ) );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( PowTest )
{
    dmemo::CommunicatorPtr comm( dmemo::Communicator::getCommunicatorPtr() );

    const IndexType n = 100;

    TestVectors vectors;

    dmemo::DistributionPtr vectorDist( new dmemo::BlockDistribution( n, comm ) );

    for ( size_t i = 0; i < vectors.size(); ++i )
    {
        Vector& v1 = *vectors[i];

        if ( ! common::isNumeric( v1.getValueType() ) )
        {
            continue;   // this test would fail for IndexType
        }

        float fillRate = 0.1f;

        v1.setSparseRandom( vectorDist, 0, fillRate, 2 );

        v1 += 2.0;   // range 2 .. 4

        VectorPtr v2Ptr( v1.copy() );
        Vector& v2 = *v2Ptr;

        v1.powExp( 2 );   // v[i] = v[i] ** 2.0
        v1.powExp( 0.5 ); 

        v1 -= v2;

        Scalar diff = v1.maxNorm();

        BOOST_CHECK( diff < Scalar( 0.0001 ) );
 
        Scalar e( common::Math::exp( 0.0 ) );

        v1.powBase( e );  // v1[i] = 2 ** v1
        v2.exp();
        v1 -= v2;

        BOOST_CHECK( diff < Scalar( 0.0001 ) );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( assign_S_VV_Test )
{
    const IndexType n = 13;

    TestVectors vectors;

    dmemo::TestDistributions dists( n );

    for ( size_t i = 0; i < vectors.size(); ++i )
    {
        VectorPtr v1 = vectors[i];

        if ( ! common::isNumeric( v1->getValueType() ) )
        {
            continue;   // this test does not work for int, uint, ....
        }

        SCAI_LOG_INFO( logger, "assign_SVV with " << *v1 )

        for ( size_t j = 0; j < dists.size(); ++j )
        {
            dmemo::DistributionPtr dist = dists[j];

            v1->allocate( dist );

            SCAI_LOG_DEBUG( logger, "dist " << j << " of " << dists.size() << ", v1 = " << *v1 )

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

    dmemo::TestDistributions dists( n );

    for ( size_t i = 0; i < vectors.size(); ++i )
    {
        VectorPtr v1 = vectors[i];

        if ( v1->getVectorKind() != Vector::DENSE )
        {
            break;
        }

        _DenseVector& dV1 = reinterpret_cast<_DenseVector&>( *v1 );

        if ( ! common::isNumeric( v1->getValueType() ) )
        {
            continue;   // this test does not work for int, uint, ....
        }

        for ( size_t j = 0; j < dists.size(); ++j )
        {
            dmemo::DistributionPtr dist = dists[j];

            dV1.setRange( dist, 3, 2 );   // only supported for dense vectors

            MatrixPtr m( Matrix::getMatrix( Matrix::CSR, v1->getValueType() ) );
            m->setIdentity( dist );

            VectorPtr v2( dV1.newVector() );

            *v2 = *m * dV1;

            // Now v1 and v2 must be equal

            *v2 -= dV1;

            BOOST_CHECK( v2->maxNorm() < Scalar( 1e-4 ) );

            *v2 = 2 * *m * dV1;

            // Now v1 and v2 must be equal

            *v2 -= 2 * dV1;

            BOOST_CHECK( v2->maxNorm() < Scalar( 1e-4 ) );

            *v2 = *m * dV1 - dV1;

            BOOST_CHECK( v2->maxNorm() < Scalar( 1e-4 ) );

            *v2 = 2 * *m * dV1 - 2 * dV1;

            BOOST_CHECK( v2->maxNorm() < Scalar( 1e-4 ) );
        }
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( assign_VM_Test )
{
    const IndexType n = 13;

    TestVectors vectors;

    dmemo::TestDistributions dists( n );

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

            MatrixPtr m( Matrix::getMatrix( Matrix::CSR, v1->getValueType() ) );
            m->setIdentity( dist );
            m->setCommunicationKind( Matrix::ASYNCHRONOUS );

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

    // we want to compare all combination of sparse/dense vectors

    TestVectors vectors1;
    TestVectors vectors2;

    dmemo::TestDistributions dists( n );

    std::srand( 1311 );   // same random numbers on all processors

    for ( size_t i = 0; i < vectors1.size(); ++i )
    {
        for ( size_t j = 0; j < vectors2.size(); ++j )
        {
            VectorPtr v1 = vectors1[i];
            VectorPtr v2 = vectors2[j];

            if ( ! common::isNumeric( v1->getValueType() ) )
            {
                continue;   // this test does not work for int, uint, ....
            }

            if ( v1->getValueType() != v2->getValueType() )
            {
                continue;   // not yet: v1 and v2 must have same type
            }

            // generate two arrays of same value type with random numbers

            common::shared_ptr<_HArray> data1( _HArray::create( v1->getValueType() ) );
            common::shared_ptr<_HArray> data2( _HArray::create( v2->getValueType() ) );

            data1->resize( n );
            data2->resize( n );

            utilskernel::HArrayUtils::setRandom( *data1, 1 );
            utilskernel::HArrayUtils::setRandom( *data2, 1 );

            v1->assign( *data1 );
            v2->assign( *data2 );

            SCAI_LOG_DEBUG( logger, "dotProduct, v1 = " << *v1 << ", v2 = " << *v2 );

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
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( scaleTest )
{
    using namespace hmemo;

    const IndexType n = 13;

    TestVectors vectors;

    dmemo::TestDistributions dists( n );

    std::srand( 1311 );   // same random numbers on all processors

    for ( size_t i = 0; i < vectors.size(); ++i )
    {
        VectorPtr v1 = vectors[i];

        if ( ! common::isNumeric( v1->getValueType() ) )
        {
            continue;   // this test does not work for int, uint, ....
        }

        common::shared_ptr<_HArray> data1( _HArray::create( v1->getValueType() ) );
        data1->resize( n );
        utilskernel::HArrayUtils::setRandom( *data1, 1 );

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
