/**
 * @file AllMatrixTest.cpp
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
 * @brief Test cases applied to each matrix class, i.e. test (virtual) methods of Matrix
 * @author Thomas Brandes
 * @date 31.08.2012
 */

#include <boost/test/unit_test.hpp>

#include <scai/lama/test/TestMacros.hpp>
#include <scai/lama/test/matrix/Matrices.hpp>

#include <scai/utilskernel/LArray.hpp>

#include <scai/hmemo/ReadAccess.hpp>
#include <scai/hmemo/HArrayRef.hpp>
#include <scai/dmemo/Distribution.hpp>

#include <scai/common/TypeTraits.hpp>
#include <scai/lama/storage/DenseStorage.hpp>

#include <scai/logging.hpp>

#include <scai/common/test/TestMacros.hpp>

using namespace scai;
using namespace lama;
using namespace dmemo;
using utilskernel::LArray;

void initMatrix( Matrix& matrix, const char* rowDistKind, const char* colDistKind )
{
    typedef SCAI_TEST_TYPE ValueType;
    const IndexType numRows = 4;
    const IndexType numColumns = 5;
    static const ValueType values[] =  { 6, 0, 7, 0, 0,
                                         0, 1, 0, 0, 0,
                                         0, 0, 9, 4, 0,
                                         2, 5, 0, 3, 8
                                       };
    hmemo::HArrayRef<ValueType> data( numRows * numColumns, values );
    DenseStorage<ValueType> denseStorage( data, numRows, numColumns );
    matrix.assign( denseStorage );
    CommunicatorPtr comm = Communicator::getCommunicatorPtr();
    // get the distributions from the factory
    DistributionPtr rowDist( Distribution::getDistributionPtr( rowDistKind, comm, numRows ) );
    DistributionPtr colDist( Distribution::getDistributionPtr( colDistKind, comm, numColumns ) );
    matrix.redistribute( rowDist, colDist );
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_SUITE( AllMatrixTest );

SCAI_LOG_DEF_LOGGER( logger, "Test.AllMatrixTest" )

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE( factoryTest )
{
    Matrices allMatrices;    // is created by factory
    size_t nFormats = Format::UNDEFINED;
    size_t nTypes   = SCAI_COMMON_COUNT_NARG( SCAI_ARITHMETIC_HOST );
    nFormats--;   // SPARSE_ASSEMBLY_STORAGE not used for a matrix
    SCAI_LOG_INFO( logger, "Test all matrices of factory to be empty, #matrices = " << allMatrices.size() )
    BOOST_CHECK_EQUAL( nTypes * nFormats, allMatrices.size() );

    for ( size_t i = 0; i < allMatrices.size(); ++i )
    {
        Matrix& matrix = *allMatrices[i];
        BOOST_CHECK_EQUAL( 0, matrix.getNumRows() );
        BOOST_CHECK_EQUAL( 0, matrix.getNumColumns() );
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE( writeAtTest )
{
    hmemo::ContextPtr context = hmemo::Context::getContextPtr();  // test context
    Matrices allMatrices( context );    // is created by factory
    SCAI_LOG_INFO( logger, "Test " << allMatrices.size() << "  matrices for typeName" )

    for ( size_t s = 0; s < allMatrices.size(); ++s )
    {
        Matrix& matrix = *allMatrices[s];
        std::ostringstream os;
        os << matrix;    // calls virtutal method writeAt for each matrix class
        BOOST_CHECK( os.str().length() > 0 );
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE( copyTest )
{
    hmemo::ContextPtr context = hmemo::Context::getContextPtr();  // test context
    // For copy we just take one arithmetic type to reduce number of test cases
    common::scalar::ScalarType stype = common::TypeTraits<SCAI_TEST_TYPE>::stype;
    Matrices allMatrices( stype, context );    // is created by factory
    SCAI_LOG_INFO( logger, "Test " << allMatrices.size() << "  matrices for copy method" )

    for ( size_t s = 0; s < allMatrices.size(); ++s )
    {
        Matrix& matrix = *allMatrices[s];
        initMatrix( matrix, "BLOCK", "BLOCK" );
        MatrixPtr copyMatrix( matrix.copy() );
        SCAI_LOG_DEBUG( logger, "copyTest: " << matrix << " with copy " << *copyMatrix );
        // verify for same matrix
        BOOST_CHECK_EQUAL( matrix.getRowDistribution(), copyMatrix->getRowDistribution() );
        BOOST_CHECK_EQUAL( matrix.getColDistribution(), copyMatrix->getColDistribution() );
        BOOST_CHECK_EQUAL( matrix.getFormat(), copyMatrix->getFormat() );
        BOOST_CHECK_EQUAL( matrix.getValueType(), copyMatrix->getValueType() );
        BOOST_CHECK_EQUAL( matrix.getContextPtr(), copyMatrix->getContextPtr() );
        Scalar maxDiff = copyMatrix->maxDiffNorm( matrix );
        // value should be exactly 0
        BOOST_CHECK_EQUAL( 0, maxDiff.getValue<float>() );
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE( l1NormTest )
{
    const IndexType N = 8;
    const Scalar scale( 2 );
    const Scalar expectedNorm = scale * Scalar( N );
    hmemo::ContextPtr context = hmemo::Context::getContextPtr();  // test context
    Matrices allMatrices( context );    // is created by factory
    SCAI_LOG_INFO( logger, "Test " << allMatrices.size() << "  matrices for l1Norm" )

    for ( size_t s = 0; s < allMatrices.size(); ++s )
    {
        Matrix& matrix = *allMatrices[s];
        matrix.setIdentity( N );
        matrix *= scale;
        Scalar l1Norm = matrix.l1Norm();
        SCAI_CHECK_CLOSE( expectedNorm, l1Norm, 0.001 );
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE( l2NormTest )
{
    const IndexType N = 8;
    const Scalar scale( 2 );
    const Scalar expectedNorm = sqrt( Scalar( N ) * scale * scale );
    hmemo::ContextPtr context = hmemo::Context::getContextPtr();  // test context
    Matrices allMatrices( context );    // is created by factory
    SCAI_LOG_INFO( logger, "Test " << allMatrices.size() << "  matrices for l2Norm" )

    for ( size_t s = 0; s < allMatrices.size(); ++s )
    {
        Matrix& matrix = *allMatrices[s];
        matrix.setIdentity( N );
        matrix *= scale;
        Scalar l2Norm = matrix.l2Norm();
        SCAI_CHECK_CLOSE( expectedNorm, l2Norm, 0.001 );
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE( maxNormTest )
{
    const IndexType N = 8;
    const Scalar scale( 2 );
    const Scalar expectedNorm = scale;
    hmemo::ContextPtr context = hmemo::Context::getContextPtr();  // test context
    Matrices allMatrices( context );    // is created by factory
    SCAI_LOG_INFO( logger, "Test " << allMatrices.size() << "  matrices for maxNorm" )

    for ( size_t s = 0; s < allMatrices.size(); ++s )
    {
        Matrix& matrix = *allMatrices[s];
        matrix.setIdentity( N );
        matrix *= scale;
        Scalar maxNorm = matrix.maxNorm();
        SCAI_CHECK_CLOSE( expectedNorm, maxNorm, 0.001 );
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE( transposeTest )
{
    hmemo::ContextPtr context = hmemo::Context::getContextPtr();  // test context
    // For transpose we just take one arithmetic type to reduce number of test cases
    // we also take same type
    common::scalar::ScalarType stype1 = common::TypeTraits<SCAI_TEST_TYPE>::stype;
    common::scalar::ScalarType stype2 = common::TypeTraits<SCAI_TEST_TYPE>::stype;
    Matrices allMatrices1( stype1, context );    // is created by factory
    Matrices allMatrices2( stype2, context );    // is created by factory
    SCAI_LOG_INFO( logger, "Test " << allMatrices1.size() << "  matrices for assignTranpose method" )

    for ( size_t s = 0; s < allMatrices1.size(); ++s )
    {
        Matrix& matrix = *allMatrices1[s];
        initMatrix( matrix, "BLOCK", "CYCLIC" );

        for ( size_t ss = 0; ss < allMatrices2.size(); ++ss )
        {
            Matrix& matrixT = *allMatrices2[ss];

            if ( matrix.getMatrixKind() != matrixT.getMatrixKind() )
            {
                continue;   // transpose for mixed sparse / dense not supported
            }

            // transpse the matrix first time
            matrixT.assignTranspose( matrix );
            SCAI_LOG_DEBUG( logger, "transposeTest: " << matrixT << " , orig is " << matrix );
            BOOST_CHECK_EQUAL( matrix.getRowDistribution(), matrixT.getColDistribution() );
            BOOST_CHECK_EQUAL( matrix.getColDistribution(), matrixT.getRowDistribution() );
            MatrixPtr matrixTT( matrix.newMatrix() );
            matrixTT->assignTranspose( matrixT );
            // verify for same matrix
            BOOST_CHECK_EQUAL( matrix.getRowDistribution(), matrixTT->getRowDistribution() );
            BOOST_CHECK_EQUAL( matrix.getColDistribution(), matrixTT->getColDistribution() );
            Scalar maxDiff = matrixTT->maxDiffNorm( matrix );

            if ( matrix.getValueType() == matrixT.getValueType() )
            {
                BOOST_CHECK_EQUAL( 0, maxDiff.getValue<SCAI_TEST_TYPE>() );
            }
            else
            {
                // we might have some conversion loss
                BOOST_CHECK( maxDiff.getValue<float>() < 0.001 );
            }
        }
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE( selfTransposeTest )
{
    hmemo::ContextPtr context = hmemo::Context::getContextPtr();  // test context
    // For transpose we just take one arithmetic type to reduce number of test cases
    common::scalar::ScalarType stype = common::TypeTraits<SCAI_TEST_TYPE>::stype;
    Matrices allMatrices( stype, context );    // is created by factory
    SCAI_LOG_INFO( logger, "Test " << allMatrices.size() << "  matrices for assignTranpose method" )

    for ( size_t s = 0; s < allMatrices.size(); ++s )
    {
        Matrix& matrix = *allMatrices[s];
        initMatrix( matrix, "BLOCK", "CYCLIC" );
        MatrixPtr copyMatrix( matrix.copy() );
        // transpse the matrix first time
        matrix.assignTranspose( matrix );
        SCAI_LOG_DEBUG( logger, "transposeTest: " << matrix << " , orig is " << *copyMatrix );
        BOOST_CHECK_EQUAL( matrix.getRowDistribution(), copyMatrix->getColDistribution() );
        BOOST_CHECK_EQUAL( matrix.getColDistribution(), copyMatrix->getRowDistribution() );
        matrix.assignTranspose( matrix );
        // verify for same matrix
        BOOST_CHECK_EQUAL( matrix.getRowDistribution(), copyMatrix->getRowDistribution() );
        BOOST_CHECK_EQUAL( matrix.getColDistribution(), copyMatrix->getColDistribution() );
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_SUITE_END();
