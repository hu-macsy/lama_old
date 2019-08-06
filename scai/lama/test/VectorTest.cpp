/**
 * @file VectorTest.cpp
 *
 * @license
 * Copyright (c) 2009-2018
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the SCAI framework LAMA.
 *
 * LAMA is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * LAMA is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
 * more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with LAMA. If not, see <http://www.gnu.org/licenses/>.
 * @endlicense
 *
 * @brief Contains generic tests for Vector objects.
 * @author Thomas Brandes
 * @date 27.07.2016
 */

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <scai/common/test/TestMacros.hpp>
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

    _TestVectors vectors;

    for ( size_t i = 0; i < vectors.size(); ++i )
    {
        _VectorPtr v = vectors[i];
        v->allocate( n );

        std::ostringstream out1;
        out1 << *v;   // same as v1->writeAt( out1 );
        BOOST_CHECK( out1.str().length() > 0 );

        std::ostringstream out2;
        v->_Vector::writeAt( out2 );
        BOOST_CHECK( out1.str().length() > 0 );

        BOOST_CHECK( out2.str() != out1.str() );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( AllocateTest, ValueType, scai_numeric_test_types )
{
    dmemo::CommunicatorPtr comm( dmemo::Communicator::getCommunicatorPtr() );

    const IndexType n = 13;

    TestVectors<ValueType> vectors;

    dmemo::DistributionPtr dist( new dmemo::BlockDistribution( n, comm ) );
    dmemo::DistributionPtr dist1( new dmemo::BlockDistribution( n + 1, comm ) );
    dmemo::DistributionPtr repDist1( new dmemo::NoDistribution( n + 1 ) );

    for ( size_t i = 0; i < vectors.size(); ++i )
    {
        VectorPtr<ValueType> v = vectors[i];

        size_t size0 = v->getMemoryUsage();

        v->allocate( n );

        size_t size1 = v->getMemoryUsage();

        if ( v->getVectorKind() == VectorKind::DENSE )
        {
            // a dense vector allocates really memory

            BOOST_CHECK( size1 >= size0 + n * common::typeSize( v->getValueType() ) );
        }
        else if ( v->getVectorKind() == VectorKind::SPARSE )
        {
            // a sparse vector does not allocate here memory

            BOOST_CHECK_EQUAL( size1, size0 );
        }

        *v = 1;
        v->redistribute( dist );
        BOOST_CHECK_EQUAL( v->getValue( n - 1 ), ValueType( 1 ) );

        BOOST_CHECK_THROW(
        {
            v->redistribute( dist1 );
        }, common::Exception );

        v->allocate( dist1 );
        *v = 2;
 

        v->redistribute( repDist1 );

        continue;

        BOOST_CHECK_EQUAL( v->getValue( n - 1 ), ValueType( 2 ) );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( SetGetTest, ValueType, scai_array_test_types )
{
    dmemo::CommunicatorPtr comm( dmemo::Communicator::getCommunicatorPtr() );

    const IndexType n = 13;

    TestVectors<ValueType> vectors;

    for ( size_t i = 0; i < vectors.size(); ++i )
    {
        Vector<ValueType>& v =*vectors[i];

        v.allocate( n );

        v = 1;

        ValueType s = v[0];
        BOOST_CHECK_EQUAL( s, ValueType( 1 ) );

        v[n-2] = 9;
        v[1] = 7;

        BOOST_CHECK_THROW(
        {
            v[n] = 1;
        }, common::Exception );

        s = v[2];
        BOOST_CHECK_EQUAL( s, ValueType( 1 ) );
        s = v[1];
        BOOST_CHECK_EQUAL( s, ValueType( 7 ) );
        v[1] = 5;
        s = v[1];
        BOOST_CHECK_EQUAL( s, ValueType( 5 ) );
        v[1] = v[2];
        s = v[1];
        BOOST_CHECK_EQUAL( s, ValueType( 1 ) );

        // Indexing of const vector has its own methods

        const Vector<ValueType>& cv = v;

        BOOST_CHECK_EQUAL( cv[1], ValueType( 1 ) );
        BOOST_CHECK_EQUAL( cv(2), ValueType( 1 ) );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( ConversionTest, ValueType, scai_numeric_test_types )
{
    typedef SCAI_TEST_TYPE OtherValueType;

    dmemo::CommunicatorPtr comm( dmemo::Communicator::getCommunicatorPtr() );

    const IndexType n = 13;

    // define some vector used for operations later

    hmemo::HArray<IndexType> indexes( { 1, 7, 11 } );
    hmemo::HArray<OtherValueType> values( { 5, 7, 9 } );

    OtherValueType zero = 1;

    SparseVector<OtherValueType> sparseVector( n, indexes, values, zero );

    auto denseVector = convert<DenseVector<OtherValueType>>( sparseVector );

    TestVectors<ValueType> vectors;

    dmemo::DistributionPtr dist = dmemo::blockDistribution( n, comm );

    RealType<ValueType> eps = common::TypeTraits<ValueType>::small();

    for ( size_t i = 0; i < vectors.size(); ++i )
    {
        Vector<ValueType>& v = *vectors[i];

        v.setSameValue( dist, 4 );

        sparseVector.redistribute( dist );
        denseVector.redistribute( dist );

        v += cast<ValueType>( denseVector );
        v -= cast<ValueType>( denseVector );
        v += cast<ValueType>( sparseVector );
        v -= cast<ValueType>( sparseVector );
        v *= cast<ValueType>( sparseVector );
        v /= cast<ValueType>( sparseVector );
        v *= cast<ValueType>( denseVector );
        v /= cast<ValueType>( denseVector );

        ValueType s = v.sum();
        ValueType expected = 4 * n;

        SCAI_LOG_INFO( logger, "sum( v ) = " << s << ", expected " << expected 
                                 << " = 4 * " << n << ", v = " << v )

        BOOST_CHECK( common::Math::abs( s - expected )  < eps );

        v = cast<ValueType>( denseVector );
        
        s = v.sum();
        expected= denseVector.sum();

        BOOST_CHECK( common::Math::abs( s - expected )  < eps );

        v = cast<ValueType>( sparseVector );
        
        s = v.sum();
        expected= sparseVector.sum();

        BOOST_CHECK( common::Math::abs( s - expected )  < eps );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( InvertTest, ValueType, scai_numeric_test_types )
{
    dmemo::CommunicatorPtr comm( dmemo::Communicator::getCommunicatorPtr() );

    const IndexType n = 13;

    TestVectors<ValueType> vectors;

    dmemo::DistributionPtr dist( new dmemo::BlockDistribution( n, comm ) );

    for ( size_t i = 0; i < vectors.size(); ++i )
    {
        Vector<ValueType>& v = *vectors[i];

        v.allocate( dist );

        v = 4;

        // v.unaryOpInPlace( common::UnaryOp::RECIPROCAL );
 
        v = 1 / v;

        ValueType s = v[n / 2 ];

        // s should be 0.25, but might not be exact

        ValueType expected = 0.25;

        RealType<ValueType> eps  = 0.00001;

        BOOST_CHECK( common::Math::abs( s - expected )  < eps );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( ConcatenateTest )
{
    typedef SCAI_TEST_TYPE ValueType;

    dmemo::CommunicatorPtr comm( dmemo::Communicator::getCommunicatorPtr() );

    const ValueType raw_values1[] = { 5, 11, 3, 2, 7, 8, 1 };
    const ValueType raw_values2[] = { 8, 14, 1, 5 };

    IndexType n1 = sizeof( raw_values1 ) / sizeof( ValueType );
    IndexType n2 = sizeof( raw_values2 ) / sizeof( ValueType );

    // build serial result by hand that is used for comparison

    DenseVector<ValueType> result;
    result.allocate( n1 + n2 );

    for ( IndexType i = 0; i < result.size(); ++i )
    {
        if ( i < n1 )
        {
            result[i] = raw_values1[i];
        }
        else
        {
            result[i] = raw_values2[i - n1];
        }
    }

    TestVectors<ValueType> vectors1;
    TestVectors<ValueType> vectors2;

    for ( size_t i1 = 0; i1 < vectors1.size(); ++i1 )
    {
        for ( size_t i2 = 0; i2 < vectors2.size(); ++i2 )
        {
            Vector<ValueType>& v1 = *vectors1[i1];
            Vector<ValueType>& v2 = *vectors2[i2];

            v1.setRawData( n1, raw_values1 );
            v2.setRawData( n2, raw_values2 );

            DenseVector<ValueType> sv;

            sv.cat( v1, v2 );

            BOOST_CHECK_EQUAL( sv.maxDiffNorm( result ), 0 );
        }
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( ConjTest, ValueType, scai_numeric_test_types )
{
    dmemo::CommunicatorPtr comm( dmemo::Communicator::getCommunicatorPtr() );

    const IndexType n = 100;

    TestVectors<ValueType> vectors;

    dmemo::DistributionPtr dist( new dmemo::BlockDistribution( n, comm ) );

    RealType<ValueType> eps = common::TypeTraits<ValueType>::small();

    for ( size_t i = 0; i < vectors.size(); ++i )
    {
        Vector<ValueType>& v = *vectors[i];

        if ( ! common::isComplex( v.getValueType() ) )
        {
            continue;   // test only useful for complex numbers
        }

        std::srand( i + 3351 );  // same values on each processor

        float fillRate = 0.2f;

        v.setSparseRandom( dist, 0, fillRate, 1 );
        
        VectorPtr<ValueType> v1Ptr( v.copy() );
        Vector<ValueType>& v1 = *v1Ptr;

        v.unaryOpInPlace( common::UnaryOp::CONJ );

        v.binaryOp( v, common::BinaryOp::MULT, v1 );  // ( a + b i ) ( a - b i ), elementwise multiplication
        
        ValueType s1 = v.sum();
        ValueType s2 = v1.dotProduct( v1 );

        SCAI_LOG_INFO( logger, "sum( v * conj(v ) = " << s1 << ", dotProduct( v, v ) = " << s2 << ", v = " << v )

        BOOST_CHECK( common::Math::abs( s1 - s2 ) < eps );

        break;
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( BinaryOpTest )
{
    typedef DefaultReal ValueType;

    const IndexType N = 3;

    auto v1 = sparseVector<ValueType>( N, 3 );
    auto v2 = sparseVector<ValueType>( N, 5 );

    v1.binaryOp( v1, common::BinaryOp::MULT, v2 );

    for ( IndexType i = 0; i < N; ++i )
    {
        ValueType s = v1[i];
        BOOST_CHECK_EQUAL( s, ValueType( 15 ) );
    }

    v1.binaryOp( v1, common::BinaryOp::SUB, 5 );

    for ( IndexType i = 0; i < N; ++i )
    {
        ValueType s = v1[i];
        BOOST_CHECK_EQUAL( s, ValueType( 10 ) );
    }
 
    v1.binaryOp( 17, common::BinaryOp::SUB, v1 );

    for ( IndexType i = 0; i < N; ++i )
    {
        ValueType s = v1[i];
        BOOST_CHECK_EQUAL( s, ValueType( 7 ) );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( BinaryOpExpTest )
{
    typedef DefaultReal ValueType;

    const IndexType N = 20;

    auto v = denseVectorLinear<ValueType>( N, -5, 1 );

    v = min( 0, v );
    v = max( v, 5 );

    for ( IndexType i = 0; i < N; ++i )
    {
        ValueType s = v[i];
        ValueType expected = common::Math::min<ValueType>( 0, i - 5 );
        expected = common::Math::max<ValueType>( expected, 5 );
        BOOST_CHECK_EQUAL( expected, s );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( ExpLogTest, ValueType, scai_numeric_test_types )
{
    dmemo::CommunicatorPtr comm( dmemo::Communicator::getCommunicatorPtr() );

    const IndexType n = 100;

    TestVectors<ValueType> vectors;

    dmemo::DistributionPtr dist( new dmemo::BlockDistribution( n, comm ) );

    for ( size_t i = 0; i < vectors.size(); ++i )
    {
        Vector<ValueType>& v1 = *vectors[i];

        float fillRate = 0.1f;

        ValueType zero = 0;

        v1.setSparseRandom( dist, zero, fillRate, 1 );

        v1 += 2;

        VectorPtr<ValueType> v2Ptr( v1.copy() );
        const Vector<ValueType>& v2 = *v2Ptr;

        v1 = exp( v1 );
        v1 = log( v1 );

        v1 -= v2;

        RealType<ValueType> diff = v1.maxNorm();
        RealType<ValueType> eps  = common::TypeTraits<ValueType>::small();

        SCAI_LOG_DEBUG( logger, "v = " << v1 << ": maxNorm( log( exp ( v ) ) - v ) = " << diff )

        BOOST_CHECK( diff < eps );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( SinCosTest, ValueType, scai_numeric_test_types )
{
    dmemo::CommunicatorPtr comm( dmemo::Communicator::getCommunicatorPtr() );

    const IndexType n = 100;

    TestVectors<ValueType> vectors;

    dmemo::DistributionPtr vectorDist( new dmemo::BlockDistribution( n, comm ) );

    for ( size_t i = 0; i < vectors.size(); ++i )
    {
        Vector<ValueType>& v1 = *vectors[i];

        float fillRate = 0.1f;

        v1.setSparseRandom( vectorDist, 0, fillRate, 1 );

        VectorPtr<ValueType> v2Ptr( v1.newVector() );
        Vector<ValueType>& v2 = *v2Ptr;

        // build:  sin(v1) * sin(v1) + cos(v2) * cos(v2) - 1, must all be 0

        v2 = cos( v1 );
        v1 = sin( v1 );

        // v1 = v1 * v1 - v2 * v2 - 1

        v1.binaryOp( v1, common::BinaryOp::MULT, v1 );
        v2.binaryOp( v2, common::BinaryOp::MULT, v2 );

        v1 += v2;
        v1 -= 1;

        RealType<ValueType> diff = v1.maxNorm();
        RealType<ValueType> eps  = 0.00001;

        BOOST_CHECK( diff < eps );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( PowTest, ValueType, scai_numeric_test_types )
{
    dmemo::CommunicatorPtr comm( dmemo::Communicator::getCommunicatorPtr() );

    const IndexType n = 100;

    TestVectors<ValueType> vectors;

    dmemo::DistributionPtr vectorDist( new dmemo::BlockDistribution( n, comm ) );

    for ( size_t i = 0; i < vectors.size(); ++i )
    {
        Vector<ValueType>& v1 = *vectors[i];

        float fillRate = 0.1f;

        ValueType zero = 0;

        IndexType bound = 2;

        v1.setSparseRandom( vectorDist, zero, fillRate, bound );

        v1 += 2.0;   // range 2 .. 4

        VectorPtr<ValueType> v2Ptr( v1.copy() );
        Vector<ValueType>& v2 = *v2Ptr;

        v1.binaryOp( v1, common::BinaryOp::POW, ValueType( 2 ) );     // v[i] = v[i] ** 2.0
        v1.binaryOp( v1, common::BinaryOp::POW, ValueType( 0.5 ) );   // v[i] = v[i] ** 0.5;

        v1 -= v2;

        RealType<ValueType> diff = v1.maxNorm();
        RealType<ValueType> eps  = 0.0001;

        BOOST_CHECK( diff < eps );
 
        ValueType e( common::Math::exp( ValueType( 0 ) ) );

        v1.binaryOp( e, common::BinaryOp::POW, v1 );  // v1[i] = e ** v1
        v2 = exp( v2 );
        v1 -= v2;

        BOOST_CHECK( diff < eps );
    }
}

/* --------------------------------------------------------------------- */

#ifdef SCAI_COMPLEX_SUPPORTED

BOOST_AUTO_TEST_CASE_TEMPLATE( ComplexTest, ValueType, scai_numeric_test_types )
{
    // skip this test if ValueType is not complex as imag would return 0

    if ( !common::isComplex( common::TypeTraits<ValueType>::stype ) )
    {
        return;
    }

    typedef RealType<ValueType> Real;
    typedef common::Complex<Real> Complex;

    dmemo::CommunicatorPtr comm( dmemo::Communicator::getCommunicatorPtr() );

    const IndexType n = 100;

    TestVectors<Complex> vectors;

    dmemo::DistributionPtr vectorDist( new dmemo::BlockDistribution( n, comm ) );

    for ( size_t i = 0; i < vectors.size(); ++i )
    {
        Vector<Complex>& complexVector = *vectors[i];

        float fillRate = 0.1f;

        Complex zero = 0;

        IndexType bound = 2;

        complexVector.setSparseRandom( vectorDist, zero, fillRate, bound );

        auto x = denseVectorEval( real ( complexVector ) );

        DenseVector<Real> y;
        y = imag( complexVector );

        auto z = denseVectorEval( complex( x, y ) );

        Real diff = complexVector.maxDiffNorm( z );

        BOOST_CHECK_EQUAL( diff, 0 );

        auto x1 = convert<DenseVector<Real>>( complexVector );

        diff = x.maxDiffNorm( x1 );

        BOOST_CHECK_EQUAL( diff, 0 );
    }
}

#endif

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( assign_S_VV_Test, ValueType, scai_numeric_test_types )
{
    const IndexType n = 13;

    TestVectors<ValueType> vectors;

    dmemo::TestDistributions dists( n );

    for ( size_t i = 0; i < vectors.size(); ++i )
    {
        Vector<ValueType>& v1 = *vectors[i];

        SCAI_LOG_INFO( logger, "assign_SVV with " << v1 )

        for ( size_t j = 0; j < dists.size(); ++j )
        {
            dmemo::DistributionPtr dist = dists[j];

            v1.allocate( dist );

            SCAI_LOG_DEBUG( logger, "dist " << j << " of " << dists.size() << ", v1 = " << v1 )

            v1 = 3;
            VectorPtr<ValueType> v2( v1.copy() );
            *v2 = 5;
            VectorPtr<ValueType> v3( v1.newVector() );
            VectorPtr<ValueType> v4( v1.newVector() );
            *v3 = v1 * *v2;
            *v4 = 3 * v1 * *v2;
            *v4 -= 2 * *v3;

            // Now v3 and v4 must be equal

            *v3 -= *v4;

            RealType<ValueType> eps = 1e-4;

            BOOST_CHECK( v3->maxNorm() < eps );
        }
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( assign_MV_Test, ValueType, scai_numeric_test_types )
{
    const IndexType n = 13;

    dmemo::TestDistributions dists( n );

    DenseVector<ValueType> dV1;

    for ( size_t j = 0; j < dists.size(); ++j )
    {
        dmemo::DistributionPtr dist = dists[j];

        IndexType bound = 5;   // random values between 0 and 5

        dV1.setRandom( dist, bound );  

        CSRSparseMatrix<ValueType> m;

        m.setIdentity( dist );

        DenseVector<ValueType> v2;
 
        SCAI_LOG_DEBUG( logger, "MV test, m = " << m << ", dV = " << dV1 << ", v2 = " << v2 )

        v2 = m * dV1;   // v1 and v2 are now equal

        v2 -= dV1;

        RealType<ValueType> eps = 1e-4;

        BOOST_CHECK( v2.maxNorm() < eps );

        v2 = 2 * m * dV1;

        // Now v1 and v2 must be equal

        v2 -= 2 * dV1;

        BOOST_CHECK( v2.maxNorm() < eps );

        v2 = m * dV1 - dV1;

        BOOST_CHECK( v2.maxNorm() < eps );

        v2 = 2 * m * dV1 - 2 * dV1;

        BOOST_CHECK( v2.maxNorm() < eps );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( assign_VM_Test, ValueType, scai_numeric_test_types )
{
    const IndexType n = 13;

    DenseVector<ValueType> v1;

    RealType<ValueType> eps = 1e-4;  // accuracy

    dmemo::TestDistributions dists( n );

    for ( size_t j = 0; j < dists.size(); ++j )
    {
        dmemo::DistributionPtr dist = dists[j];

        v1.allocate( dist );

        v1 = ValueType( 3 );

        CSRSparseMatrix<ValueType> m;

        m.setIdentity( dist );
        m.setCommunicationKind( SyncKind::ASYNC_LOCAL );

        VectorPtr<ValueType> v2( v1.newVector() );

        *v2 = transpose( m ) * v1;

        // Now v1 and v2 must be equal

        *v2 -= v1;

        BOOST_CHECK( v2->maxNorm() < eps );

        *v2 = 2 * transpose( m ) * v1;

        // Now 2 * v1 and v2 must be equal

        *v2 -= 2 * v1;

        BOOST_CHECK( v2->maxNorm() < eps );

        *v2 = transpose( m ) * v1 - v1;

        BOOST_CHECK( v2->maxNorm() < eps );

        *v2 = 2 * transpose( m ) * v1 - 2 * v1;

        BOOST_CHECK( v2->maxNorm() < eps );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( dotProductTest, ValueType, scai_numeric_test_types )
{
    // This test does not verify the correctness of the dot product
    // but that it is the same, either computed replicated or distributed.

    using namespace hmemo;

    const IndexType n = 13;

    // we want to compare all combination of sparse/dense vectors

    TestVectors<ValueType> vectors1;
    TestVectors<ValueType> vectors2;

    dmemo::TestDistributions dists( n );

    std::srand( 1311 );   // same random numbers on all processors

    for ( size_t i = 0; i < vectors1.size(); ++i )
    {
        for ( size_t j = 0; j < vectors2.size(); ++j )
        {
            Vector<ValueType>& v1 = *vectors1[i];
            Vector<ValueType>& v2 = *vectors2[j];

            // generate two arrays of same value type with random numbers

            HArray<ValueType> data1( n );
            HArray<ValueType> data2( n );

            utilskernel::HArrayUtils::setRandom( data1, 1 );
            utilskernel::HArrayUtils::setRandom( data2, 1 );

            v1.assign( data1 );
            v2.assign( data2 );

            SCAI_LOG_DEBUG( logger, "dotProduct, v1 = " << v1 << ", v2 = " << v2 );

            ValueType dotp = v1.dotProduct( v2 );  // replicated computation

            for ( size_t j = 0; j < dists.size(); ++j )
            {
                // now compute the dot product with distributed vectors

                dmemo::DistributionPtr dist = dists[j];

                v1.redistribute( dist );
                v2.redistribute( dist );

                ValueType distDotp = v1.dotProduct( v2 );

                // we cannot check for equality due to different rounding errors

                RealType<ValueType> eps = 0.0001;

                BOOST_CHECK( common::Math::abs( dotp - distDotp ) < eps );
            }
        }
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( scaleTest, ValueType, scai_numeric_test_types )
{
    using namespace hmemo;

    const IndexType n = 13;

    TestVectors<ValueType> vectors;

    dmemo::TestDistributions dists( n );

    std::srand( 1311 );   // same random numbers on all processors

    DenseVector<ValueType> v1;

    HArray<ValueType> data1( n );
    utilskernel::HArrayUtils::setRandom( data1, 1 );

    v1.assign( data1 );

    for ( size_t j = 0; j < dists.size(); ++j )
    {
        v1.redistribute( dists[j] );

        CSRSparseMatrix<ValueType> m;

        m.setIdentity( dists[j] );
        m.setDiagonal( v1 );

        DenseVector<ValueType> v2( v1 );

        v2 = v1;

        v2 += v1;
        v2 *= 2;
        v2 /= 2;
        v2 -= v1;

        v2 = v1 * v2;   // is v2 * v2 elementwise

        DenseVector<ValueType> v3;

        v3 = m * v1 - v2;   // is v1 * v1

        RealType<ValueType> eps = 0.0001;
        BOOST_CHECK( v3.maxNorm() < eps );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( allTest, ValueType, scai_numeric_test_types )
{
    using namespace hmemo;

    const IndexType n = 5;

    const DefaultReal denseData[] = { 0, 3, 3, 5, 2 };

    const DefaultReal zero = 3;
    const IndexType sparseIndexes[] = { 0, 3, 4 };
    const DefaultReal sparseData[] = { 0, 5, 2 };
    const IndexType nnz = 3;

    // we want to compare all combination of sparse/dense vectors

    TestVectors<ValueType> vectors1;
    TestVectors<ValueType> vectors2;

    RealType<ValueType> eps = common::TypeTraits<ValueType>::small();

    for ( size_t i = 0; i < vectors1.size(); ++i )
    {
        Vector<ValueType>& v1 = *vectors1[i];

        v1.setSparseRawData( n, 3, sparseIndexes, sparseData, ValueType( 0 ) );

        BOOST_CHECK( v1.all( common::CompareOp::GE, 0 ) );

        v1.setRawData( n, denseData );

        BOOST_CHECK( v1.all( common::CompareOp::GE, 0 ) );

        for ( size_t j = 0; j < vectors2.size(); ++j )
        {
            Vector<ValueType>& v2 = *vectors2[j];

            v2.setSparseRawData( n, nnz, sparseIndexes, sparseData, zero );

            BOOST_CHECK( v1.maxDiffNorm( v2 ) < eps );
        }
    }
}

BOOST_AUTO_TEST_SUITE_END();
