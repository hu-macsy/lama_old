/**
 * @file AllMatrixTest.cpp
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
 * @brief Test cases applied to each matrix class, i.e. test (virtual) methods of _Matrix
 * @author Thomas Brandes
 * @date 31.08.2012
 */

#include <boost/test/unit_test.hpp>

#include <scai/common/test/TestMacros.hpp>
#include <scai/lama/test/matrix/Matrices.hpp>


#include <scai/hmemo/ReadAccess.hpp>
#include <scai/hmemo/HArrayRef.hpp>

#include <scai/dmemo/Distribution.hpp>
#include <scai/dmemo/BlockDistribution.hpp>
#include <scai/dmemo/RedistributePlan.hpp>
#include <scai/dmemo/test/TestDistributions.hpp>

#include <scai/common/TypeTraits.hpp>
#include <scai/common/Math.hpp>

#include <scai/lama/storage/DenseStorage.hpp>
#include <scai/lama/Scalar.hpp>
#include <scai/lama/matrix/CSRSparseMatrix.hpp>
#include <scai/lama/matrix/DenseMatrix.hpp>
#include <scai/lama/matutils/MatrixCreator.hpp>
#include <scai/sparsekernel/COOUtils.hpp>

#include <scai/logging.hpp>


using namespace scai;
using namespace lama;
using namespace dmemo;
using hmemo::HArray;

using boost::test_tools::per_element;

void initMatrix( _Matrix& matrix, const char* rowDistKind, const char* colDistKind )
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
    DenseStorage<ValueType> denseStorage( numRows, numColumns, data );
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
    _Matrices allMatrices;    // is created by factory
    size_t nFormats = static_cast<size_t>( Format::UNDEFINED );
    size_t nTypes   = SCAI_COMMON_COUNT_NARG( SCAI_NUMERIC_TYPES_HOST );
    nFormats--;   // STENCIL_STORAGE not used for a matrix in factory
    nFormats--;   // AssemblyStorage not used for a matrix
    SCAI_LOG_INFO( logger, "Test all matrices of factory to be empty, #matrices = " << allMatrices.size() )
    BOOST_CHECK_EQUAL( nTypes * nFormats, allMatrices.size() );

    for ( size_t i = 0; i < allMatrices.size(); ++i )
    {
        _Matrix& matrix = *allMatrices[i];
        BOOST_CHECK_EQUAL( IndexType( 0 ), matrix.getNumRows() );
        BOOST_CHECK_EQUAL( IndexType( 0 ), matrix.getNumColumns() );
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE( writeAtTest )
{
    hmemo::ContextPtr context = hmemo::Context::getContextPtr();  // test context
    _Matrices allMatrices( context );    // is created by factory
    SCAI_LOG_INFO( logger, "Test " << allMatrices.size() << "  matrices for typeName" )

    for ( size_t s = 0; s < allMatrices.size(); ++s )
    {
        _Matrix& matrix = *allMatrices[s];

        std::ostringstream os;
        os << matrix;    // calls virtutal method writeAt for each matrix class
        BOOST_CHECK( os.str().length() > 0 );

        // print the different enum value, length() > 1 makes sure that not an int is printed

        std::ostringstream os1;
        os1 << matrix.getMatrixKind();
        BOOST_CHECK( os1.str().length() > 1 );

        std::ostringstream os2;
        os2 << matrix.getCommunicationKind();
        BOOST_CHECK( os2.str().length() > 1 );

        std::ostringstream os3;
        os3 << matrix.getFormat();
        BOOST_CHECK( os3.str().length() > 1 );

        std::ostringstream os4;
        os4 << matrix.getTypeName();
        BOOST_CHECK( os4.str().length() > 1 );
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE( setContextTest )
{
    hmemo::ContextPtr context = hmemo::Context::getContextPtr();  // test context

    hmemo::ContextPtr nullContext = hmemo::ContextPtr();   // null context

    _Matrices allMatrices( context );    // is created by factory

    SCAI_LOG_INFO( logger, "Test " << allMatrices.size() << "  matrices for setContext" )

    for ( size_t s = 0; s < allMatrices.size(); ++s )
    {
        _Matrix& matrix = *allMatrices[s];

        BOOST_CHECK_THROW(
        {
            matrix.setContextPtr( nullContext );
        }, common::Exception );

        matrix.setContextPtr( context );

        BOOST_CHECK_EQUAL( context.get(), matrix.getContextPtr().get() );
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE_TEMPLATE( copyTest, ValueType, scai_numeric_test_types )
{
    hmemo::ContextPtr context = hmemo::Context::getContextPtr();  // test context

    Matrices<ValueType> testMatrices( context );    // is created by factory

    SCAI_LOG_INFO( logger, "Test " << testMatrices.size() << "  matrices for copy method" )

    for ( size_t s = 0; s < testMatrices.size(); ++s )
    {
        Matrix<ValueType>& matrix = *testMatrices[s];
        initMatrix( matrix, "BLOCK", "BLOCK" );
        MatrixPtr<ValueType> copyMatrix( matrix.copy() );
        // verify for same matrix
        BOOST_CHECK_EQUAL( matrix.getRowDistribution(), copyMatrix->getRowDistribution() );
        BOOST_CHECK_EQUAL( matrix.getColDistribution(), copyMatrix->getColDistribution() );
        BOOST_CHECK_EQUAL( matrix.getFormat(), copyMatrix->getFormat() );
        BOOST_CHECK_EQUAL( matrix.getValueType(), copyMatrix->getValueType() );
        BOOST_CHECK_EQUAL( matrix.getContextPtr(), copyMatrix->getContextPtr() );
        ValueType maxDiff = copyMatrix->maxDiffNorm( matrix );
        // value should be exactly 0
        BOOST_CHECK_EQUAL( ValueType( 0 ), maxDiff );
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE_TEMPLATE( l1NormTest, ValueType, scai_numeric_test_types )
{
    const IndexType N = 8;
    const ValueType scale = 2;
    const RealType<ValueType> expectedNorm = scale * ValueType( N );
    const RealType<ValueType> eps = common::TypeTraits<ValueType>::small();

    hmemo::ContextPtr context = hmemo::Context::getContextPtr();  // test context
    Matrices<ValueType> allMatrices( context );    // is created by factory
    SCAI_LOG_INFO( logger, "Test " << allMatrices.size() << "  matrices for l1Norm" )

    for ( size_t s = 0; s < allMatrices.size(); ++s )
    {
        Matrix<ValueType>& matrix = *allMatrices[s];
        matrix.setIdentity( N );
        matrix *= scale;
        SCAI_LOG_DEBUG( logger, "Test l1Norm for this matrix: " << matrix )
        RealType<ValueType> l1Norm = matrix.l1Norm();
        BOOST_CHECK( common::Math::abs( expectedNorm - l1Norm ) <  eps );
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE_TEMPLATE( l2NormTest, ValueType, scai_numeric_test_types )
{
    const IndexType N = 8;
    const ValueType scale = 2;
    const RealType<ValueType> expectedNorm = common::Math::sqrt( ValueType( N ) * scale * scale );

    hmemo::ContextPtr context = hmemo::Context::getContextPtr();  // test context

    Matrices<ValueType> allMatrices( context );    // is created by factory

    SCAI_LOG_DEBUG( logger, "Test " << allMatrices.size() << "  matrices for l2Norm" )

    for ( size_t s = 0; s < allMatrices.size(); ++s )
    {
        Matrix<ValueType>& matrix = *allMatrices[s];
        matrix.setIdentity( N );
        matrix *= scale;
        SCAI_LOG_DEBUG( logger, "Test l2Norm for this matrix: " << matrix )
        RealType<ValueType> l2Norm = matrix.l2Norm();
        BOOST_CHECK_CLOSE( expectedNorm, l2Norm, 0.001 );
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE_TEMPLATE( maxNormTest, ValueType, scai_numeric_test_types )
{
    const IndexType N = 8;
    const ValueType scale =  2;
    const RealType<ValueType> expectedNorm = scale;
    hmemo::ContextPtr context = hmemo::Context::getContextPtr();  // test context
    Matrices<ValueType> allMatrices( context );    // is created by factory
    SCAI_LOG_INFO( logger, "Test " << allMatrices.size() << "  matrices for maxNorm" )

    for ( size_t s = 0; s < allMatrices.size(); ++s )
    {
        Matrix<ValueType>& matrix = *allMatrices[s];
        matrix.setIdentity( N );
        matrix *= scale;
        SCAI_LOG_DEBUG( logger, "Test maxNorm for this matrix: " << matrix )
        RealType<ValueType> maxNorm = matrix.maxNorm();
        BOOST_CHECK_CLOSE( expectedNorm, maxNorm, 0.01 );
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE_TEMPLATE( scaleTest, ValueType, scai_numeric_test_types )
{
    const IndexType M = 10;
    const IndexType N = 8;

    hmemo::ContextPtr context = hmemo::Context::getContextPtr();  // test context

    Matrices<ValueType> allMatrices( context );    // is created by factory

    SCAI_LOG_INFO( logger, "Test " << allMatrices.size() << "  matrices for scale" )

    std::srand( 10113 );  // same random numbers on each processor

    auto input = zero<CSRSparseMatrix<ValueType>>( M, N );

    MatrixCreator::fillRandom( input, 0.5f );

    auto scaleY = denseVectorLinear<ValueType>( M, 1, 1 );

    CSRSparseMatrix<ValueType> output( input );
    output.scaleRows( scaleY );
    output.scale( 0.5 );

    TestDistributions rowDist( M );
    TestDistributions colDist( N );

    for ( size_t s = 0; s < allMatrices.size(); ++s )
    {
        Matrix<ValueType>& matrix = *allMatrices[s];

        for ( size_t i = 0; i < rowDist.size(); ++i )
        {
            scaleY.redistribute( rowDist[i] );

            for ( size_t j = 0; j < colDist.size(); ++j )
            {
                matrix = input;   // serial matrix

                matrix.redistribute( rowDist[i], colDist[j] );

                // now scale parallel

                matrix.scaleRows( scaleY );
                matrix.scale( 0.5 );

                // verify correct results by replication of matrix

                matrix.redistribute( input.getRowDistributionPtr(), input.getColDistributionPtr() );

                RealType<ValueType> diff = output.maxDiffNorm( matrix );

                SCAI_LOG_DEBUG( logger, "diff = " << diff << ", matrix = " << matrix )

                // there should be no rounding errors, so we check here for exact results

                BOOST_CHECK_EQUAL( diff, 0 );
            }
        }
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE( scaleColumnsTest )
{
    typedef DefaultReal ValueType;    // as we test here commmunication patterns, only one ValueType needed

    const IndexType M = 10;
    const IndexType N = 8;

    hmemo::ContextPtr context = hmemo::Context::getContextPtr();  // test context

    Matrices<ValueType> allMatrices( context );    // is created by factory

    SCAI_LOG_INFO( logger, "Test " << allMatrices.size() << "  matrices for scale" )

    std::srand( 10113 );  // same random numbers on each processor

    auto input = zero<CSRSparseMatrix<ValueType>>( M, N );

    MatrixCreator::fillRandom( input, 0.5f );

    auto scaleY = denseVectorLinear<ValueType>( N, 1, 1 );

    CSRSparseMatrix<ValueType> output( input );
    output.scaleColumns( scaleY );

    TestDistributions rowDist( M );
    TestDistributions colDist( N );

    for ( size_t s = 0; s < allMatrices.size(); ++s )
    {
        Matrix<ValueType>& matrix = *allMatrices[s];

        for ( size_t i = 0; i < rowDist.size(); ++i )
        {
            for ( size_t j = 0; j < colDist.size(); ++j )
            {
                matrix = input;   // serial matrix

                matrix.redistribute( rowDist[i], colDist[j] );
                scaleY.redistribute( colDist[j] );

                // now scale with distributed vector

                matrix.scaleColumns( scaleY );

                // verify correct results by replication of matrix

                matrix.redistribute( input.getRowDistributionPtr(), input.getColDistributionPtr() );

                RealType<ValueType> diff = output.maxDiffNorm( matrix );

                SCAI_LOG_DEBUG( logger, "diff = " << diff << ", matrix = " << matrix )

                // there should be no rounding errors, so we check here for exact results

                BOOST_CHECK_EQUAL( diff, 0 );
            }
        }
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE( transposeTest )
{
    hmemo::ContextPtr context = hmemo::Context::getContextPtr();  // test context

    // For transpose we just take one arithmetic type to reduce number of test cases
    // we also take same type

    typedef SCAI_TEST_TYPE ValueType;

    Matrices<ValueType> allMatrices1( context ); // one instance of each supported matrix of ValueType
    Matrices<ValueType> allMatrices2( context ); // one instance of each supported matrix of ValueType

    SCAI_LOG_INFO( logger, "Test " << allMatrices1.size() << "  matrices for assignTranpose method" )

    for ( size_t s = 0; s < allMatrices1.size(); ++s )
    {
        Matrix<ValueType>& matrix = *allMatrices1[s];

        initMatrix( matrix, "BLOCK", "CYCLIC" );

        for ( size_t ss = 0; ss < allMatrices2.size(); ++ss )
        {
            Matrix<ValueType>& matrixT = *allMatrices2[ss];

            if ( matrix.getMatrixKind() != matrixT.getMatrixKind() )
            {
                continue;   // transpose for mixed sparse / dense not supported
            }

            // transpse the matrix first time

            SCAI_LOG_DEBUG( logger, "transposeTest: orig is " << matrix );
            matrixT.assignTranspose( matrix );
            BOOST_CHECK_EQUAL( matrix.getRowDistribution(), matrixT.getColDistribution() );
            BOOST_CHECK_EQUAL( matrix.getColDistribution(), matrixT.getRowDistribution() );
            MatrixPtr<ValueType> matrixTT( matrix.newMatrix() );
            matrixTT->assignTranspose( matrixT );
            // verify for same matrix
            BOOST_CHECK_EQUAL( matrix.getRowDistribution(), matrixTT->getRowDistribution() );
            BOOST_CHECK_EQUAL( matrix.getColDistribution(), matrixTT->getColDistribution() );
            ValueType maxDiff = matrixTT->maxDiffNorm( matrix );

            BOOST_CHECK_EQUAL( 0, maxDiff );
        }
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE( selfTransposeTest )
{
    hmemo::ContextPtr context = hmemo::Context::getContextPtr();  // test context

    typedef SCAI_TEST_TYPE ValueType;

    // For transpose we just take one arithmetic type to reduce number of test cases

    Matrices<ValueType> allMatrices( context );    // is created by factory
    SCAI_LOG_INFO( logger, "Test " << allMatrices.size() << "  matrices for assignTranpose method" )

    for ( size_t s = 0; s < allMatrices.size(); ++s )
    {
        Matrix<ValueType>& matrix = *allMatrices[s];
        initMatrix( matrix, "BLOCK", "CYCLIC" );
        MatrixPtr<ValueType> copyMatrix( matrix.copy() );
        // transpse the matrix first time
        matrix = transpose( matrix );
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

BOOST_AUTO_TEST_CASE_TEMPLATE( assignAddTest, ValueType, scai_numeric_test_types )
{
    // ToDo: test fails due to alias for ELLStorage;

    hmemo::ContextPtr context = hmemo::Context::getContextPtr();  // test context

    Matrices<ValueType> testMatrices( context );

    SCAI_LOG_INFO( logger, "Test " << testMatrices.size() << "  matrices for assign operator tests" )

    for ( size_t s = 0; s < testMatrices.size(); ++s )
    {
        Matrix<ValueType>& matrix1 = *testMatrices[s];

        initMatrix( matrix1, "BLOCK", "NO" );

        MatrixPtr<ValueType> matrix2Ptr( matrix1.copy() );
        MatrixPtr<ValueType> matrix3Ptr( matrix1.newMatrix() );

        Matrix<ValueType>& matrix2 = *matrix2Ptr;
        Matrix<ValueType>& matrix3 = *matrix3Ptr;

        matrix3 = matrix1 + matrix2;

        BOOST_CHECK_EQUAL( matrix3.getRowDistribution(), matrix1.getRowDistribution() );
        BOOST_CHECK_EQUAL( matrix3.getColDistribution(), matrix1.getColDistribution() );

        matrix3 = matrix3 - matrix2;
        matrix1 += 2 * matrix3;
        matrix1 -= matrix3;
        matrix1 += matrix3;
        matrix1 -= 3 * matrix3;

        BOOST_CHECK_EQUAL( ValueType( 0 ), matrix1.maxNorm() );

        matrix1 = matrix2;

        BOOST_CHECK_EQUAL( matrix1.getRowDistributionPtr(), matrix2.getRowDistributionPtr() );
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE_TEMPLATE( assignMultTest, ValueType, scai_numeric_test_types )
{
    hmemo::ContextPtr context = hmemo::Context::getContextPtr();  // test context

    Matrices<ValueType> allMatrices( context );    // is created by factory

    SCAI_LOG_INFO( logger, "Test " << allMatrices.size() << "  matrices for assign operator tests" )

    for ( size_t s = 0; s < allMatrices.size(); ++s )
    {
        Matrix<ValueType>& matrix1 = *allMatrices[s];

        if ( matrix1.getMatrixKind() == MatrixKind::DENSE )
        {
            continue;
        }

        initMatrix( matrix1, "BLOCK", "NO" );

        CSRSparseMatrix<ValueType> unityLeft;
        CSRSparseMatrix<ValueType> unityRight;

        unityLeft.setIdentity( matrix1.getRowDistributionPtr() );
        unityRight.setIdentity( matrix1.getColDistributionPtr() );

        MatrixPtr<ValueType> matrix2Ptr( matrix1.newMatrix() );

        Matrix<ValueType>& matrix2 = *matrix2Ptr;

        matrix2 = matrix1 * unityRight;   // not for Dense

        BOOST_CHECK_EQUAL( matrix2.getRowDistributionPtr(), matrix1.getRowDistributionPtr() );

        matrix2 = matrix1 * unityRight + matrix1;

        BOOST_CHECK_EQUAL( matrix2.getRowDistributionPtr(), matrix1.getRowDistributionPtr() );
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE_TEMPLATE( checkSymmetryTest, ValueType, scai_numeric_test_types )
{
    hmemo::ContextPtr context = hmemo::Context::getContextPtr();  // test context

    Matrices<ValueType> allMatrices( context );    // is created by factory

    SCAI_LOG_INFO( logger, "Test " << allMatrices.size() << "  matrices for checkSymmetry" )

    for ( size_t s = 0; s < allMatrices.size(); ++s )
    {
        Matrix<ValueType>& matrix = *allMatrices[s];

        matrix.setIdentity( 5 );
        BOOST_CHECK( matrix.checkSymmetry() );

        MatrixCreator::buildPoisson2D( matrix, 5, 3, 3 );
        BOOST_CHECK( matrix.checkSymmetry() );
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE_TEMPLATE( diagonalTest, ValueType, scai_numeric_test_types )
{
    const IndexType n1 = 3;
    const IndexType n2 = 4;

    hmemo::ContextPtr context = hmemo::Context::getContextPtr();  // test context

    Matrices<ValueType> allMatrices( context );    // is created by factory

    TestDistributions testDistributions( n1 * n2 );

    SCAI_LOG_INFO( logger, "Test " << allMatrices.size() << "  matrices for checkSymmetry" )

    for ( size_t s = 0; s < allMatrices.size(); ++s )
    {
        Matrix<ValueType>& matrix = *allMatrices[s];

        matrix.setCommunicationKind( SyncKind::SYNCHRONOUS );

        RealType<ValueType> eps = 0.0001;

        for ( size_t i = 0; i < testDistributions.size(); ++i )
        {
            DistributionPtr dist = testDistributions[i];

            matrix.clear();

            MatrixCreator::buildPoisson2D( matrix, 5, n1, n2 );

            matrix.redistribute( dist, dist );

            SCAI_LOG_DEBUG( logger, "diagonalTest for " << matrix )

            DenseVector<ValueType> x;
            DenseVector<ValueType> y1;
            DenseVector<ValueType> y2;
            DenseVector<ValueType> d;

            x.allocate( matrix.getColDistributionPtr() );
            x  = ValueType( 1 );

            y1 = matrix * x;

            matrix.getDiagonal( d );

            BOOST_CHECK_EQUAL( d.getDistribution(), matrix.getRowDistribution() );

            matrix.setDiagonal( 0 );

            // Now we can prove y2 = matrix * x + diagonal must be same as y1

            y2 = matrix * x + d;
            y2 = y2 - y1;

            BOOST_CHECK( y2.maxNorm() < eps );

            // Write back modified diagonal and check result

            d += 1;

            matrix.setDiagonal( d );

            y2 = matrix * x - y1;
            y2 += -1;

            BOOST_CHECK( y2.maxNorm() < eps );
        }
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE_TEMPLATE( getRowTest, ValueType, scai_numeric_test_types )
{
    const IndexType nRows = 10;
    const IndexType nCols = 8;

    hmemo::ContextPtr context = hmemo::Context::getContextPtr();  // test context

    auto csr = zero<CSRSparseMatrix<DefaultReal>>( nRows, nCols );
    MatrixCreator::fillRandom( csr, 0.1f );

    TestDistributions rowDistributions( nRows );
    TestDistributions colDistributions( nCols );

    RealType<ValueType> eps = common::TypeTraits<ValueType>::small();

    for ( size_t i = 0; i < rowDistributions.size(); ++i )
    {
        DistributionPtr rowDist = rowDistributions[i];

        for ( size_t j = 0; j < colDistributions.size(); ++j )
        {
            DistributionPtr colDist = colDistributions[j];

            Matrices<ValueType> allMatrices( context );    // is created by factory

            for ( size_t s = 0; s < allMatrices.size(); ++s )
            {
                Matrix<ValueType>& matrix = *allMatrices[s];

                matrix = cast<ValueType>( csr );  // convert DefaultReal->ValueType

                matrix.redistribute( rowDist, colDist );

                SCAI_LOG_INFO( logger, "getRowTest for this matrix: " << matrix )

                // get each row and subtract it

                DenseVector<ValueType> row;

                for ( IndexType iRow = 0; iRow < matrix.getNumRows(); ++iRow )
                {
                    matrix.getRow( row, iRow );
                    // BOOST_REQUIRE_EQUAL( nCols, row->size() );
                    // BOOST_REQUIRE( row->isConsistent() );
                    matrix.setRow( row, iRow, common::BinaryOp::SUB );
                }

                // the final matrix should be zero

                BOOST_CHECK( matrix.maxNorm() < eps );
            }
        }
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE_TEMPLATE( reduceTest, ValueType, scai_numeric_test_types )
{
    const IndexType nRows = 10;
    const IndexType nCols = 8;

    hmemo::ContextPtr context = hmemo::Context::getContextPtr();  // test context

    auto csr = zero<CSRSparseMatrix<ValueType>>( nRows, nCols );
    MatrixCreator::fillRandom( csr, 0.1f );

    common::BinaryOp reduceOp = common::BinaryOp::ADD;
    common::UnaryOp   elemOp   = common::UnaryOp::SQR;

    TestDistributions rowDistributions( nRows );
    TestDistributions colDistributions( nCols );

    RealType<ValueType> eps = common::TypeTraits<ValueType>::small();

    for ( IndexType dim = 0; dim < 2; ++dim )
    {
        DenseVector<ValueType> sRow;
        csr.reduce( sRow, dim, reduceOp, elemOp );

        for ( size_t i = 0; i < rowDistributions.size(); ++i )
        {
            DistributionPtr rowDist = rowDistributions[i];

            for ( size_t j = 0; j < colDistributions.size(); ++j )
            {
                DistributionPtr colDist = colDistributions[j];

                Matrices<ValueType> allMatrices( context );    // is created by factory

                for ( size_t s = 0; s < allMatrices.size(); ++s )
                {
                    Matrix<ValueType>& matrix = *allMatrices[s];

                    matrix = csr;   // format + type conversion

                    matrix.redistribute( rowDist, colDist );

                    DenseVector<ValueType> row;

                    // reduce on the parallel matrix

                    matrix.reduce( row, dim, reduceOp, elemOp );

                    if ( dim == 0 )
                    {
                        BOOST_CHECK_EQUAL( row.size(), nRows );
                    }
                    else
                    {
                        BOOST_CHECK_EQUAL( row.size(), nCols );
                    }

                    row.redistribute( sRow.getDistributionPtr() );

                    if ( ! ( row.maxDiffNorm( sRow ) < eps ) )
                    {
                        SCAI_LOG_ERROR( logger,  "i = " << i << ", j = " << j << ", s = " << s
                                        << ", fail for this matrix: " << matrix << std::endl );
                    }

                    BOOST_CHECK( row.maxDiffNorm( sRow ) < eps );
                }
            }
        }
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE_TEMPLATE( getColTest, ValueType, scai_numeric_test_types )
{
    const IndexType nRows = 10;
    const IndexType nCols = 8;

    hmemo::ContextPtr context = hmemo::Context::getContextPtr();  // test context

    auto csr = zero<CSRSparseMatrix<DefaultReal>>( nRows, nCols );
    MatrixCreator::fillRandom( csr, 0.1f );

    TestDistributions rowDistributions( nRows );
    TestDistributions colDistributions( nCols );

    RealType<ValueType> eps = common::TypeTraits<ValueType>::small();

    for ( size_t i = 0; i < rowDistributions.size(); ++i )
    {
        DistributionPtr rowDist = rowDistributions[i];

        for ( size_t j = 0; j < colDistributions.size(); ++j )
        {
            DistributionPtr colDist = colDistributions[j];

            Matrices<ValueType> allMatrices( context );    // is created by factory

            for ( size_t s = 0; s < allMatrices.size(); ++s )
            {
                Matrix<ValueType>& matrix = *allMatrices[s];

                matrix = cast<ValueType>( csr );

                matrix.redistribute( rowDist, colDist );

                // get each row and subtract it

                DenseVector<ValueType> col;

                SCAI_LOG_INFO( logger, "getColTest for this matrix: " << matrix )

                for ( IndexType iCol = 0; iCol < matrix.getNumColumns(); ++iCol )
                {
                    matrix.getColumn( col, iCol );
                    matrix.setColumn( col, iCol, common::BinaryOp::SUB );
                }

                // the final matrix should be zero

                BOOST_CHECK( matrix.maxNorm() < eps );
            }
        }
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE_TEMPLATE( getTest, ValueType, scai_numeric_test_types )
{
    const IndexType nRows = 10;
    const IndexType nCols = 8;

    hmemo::ContextPtr context = hmemo::Context::getContextPtr();  // test context

    auto csr = zero<CSRSparseMatrix<DefaultReal>>( nRows, nCols );

    MatrixCreator::fillRandom( csr, 0.1f );

    TestDistributions rowDistributions( nRows );
    TestDistributions colDistributions( nCols );

    RealType<ValueType> eps = 0.0001;

    for ( size_t rd = 0; rd < rowDistributions.size(); ++rd )
    {
        DistributionPtr rowDist = rowDistributions[rd];

        for ( size_t cd = 0; cd < colDistributions.size(); ++cd )
        {
            DistributionPtr colDist = colDistributions[cd];

            Matrices<ValueType> allMatrices( context );    // is created by factory

            for ( size_t s = 0; s < allMatrices.size(); ++s )
            {
                Matrix<ValueType>& matrix = *allMatrices[s];

                matrix = cast<ValueType>( csr );

                matrix.redistribute( rowDist, colDist );

                SCAI_LOG_DEBUG( logger, "getTest for this matrix: " << matrix << ", max = " << matrix.maxNorm() )

                // get each element and subract it

                for ( IndexType iRow = 0; iRow < matrix.getNumRows(); ++iRow )
                {
                    for ( IndexType jCol = 0; jCol < matrix.getNumColumns(); ++jCol )
                    {
                        ValueType s1 = matrix.getValue( iRow, jCol );
                        ValueType s2 = static_cast<ValueType>( csr.getValue( iRow, jCol ) );

                        RealType<ValueType> diff = common::Math::abs( s1 - s2 );

                        BOOST_CHECK( diff < eps );
                    }
                }
            }
        }
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE_TEMPLATE( getSetTest, ValueType, scai_numeric_test_types )
{
    const IndexType nRows = 10;
    const IndexType nCols = 8;

    hmemo::ContextPtr context = hmemo::Context::getContextPtr();  // test context

    auto csr = zero<CSRSparseMatrix<ValueType>>( nRows, nCols );

    MatrixCreator::fillRandom( csr, 0.1f );

    TestDistributions rowDistributions( nRows );
    TestDistributions colDistributions( nCols );

    RealType<ValueType> eps = common::TypeTraits<ValueType>::small();

    for ( size_t rd = 0; rd < rowDistributions.size(); ++rd )
    {
        DistributionPtr rowDist = rowDistributions[rd];

        for ( size_t cd = 0; cd < colDistributions.size(); ++cd )
        {
            DistributionPtr colDist = colDistributions[cd];

            Matrices<ValueType> allMatrices( context );    // is created by factory

            for ( size_t s = 0; s < allMatrices.size(); ++s )
            {
                Matrix<ValueType>& matrix = *allMatrices[s];

                matrix = csr;

                matrix.redistribute( rowDist, colDist );

                SCAI_LOG_DEBUG( logger, "getSetTest for this matrix: " << matrix << ", max = " << matrix.maxNorm() )

                // get each element and subract it

                for ( IndexType iRow = 0; iRow < matrix.getNumRows(); ++iRow )
                {
                    for ( IndexType jCol = 0; jCol < matrix.getNumColumns(); ++jCol )
                    {
                        ValueType s = matrix.getValue( iRow, jCol );

                        if ( s != ValueType( 0 ) )
                        {
                            matrix.setValue( iRow, jCol, s, common::BinaryOp::SUB );
                        }
                    }
                }

                // the final matrix should be zero

                BOOST_CHECK( matrix.maxNorm() < eps );
            }
        }
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE_TEMPLATE( redistributeTest, ValueType, scai_numeric_test_types )
{
    const IndexType n = 15;   // only square matrices here

    hmemo::ContextPtr context = hmemo::Context::getContextPtr();  // test context

    auto csr = zero<CSRSparseMatrix<DefaultReal>>( n, n );

    MatrixCreator::fillRandom( csr, 0.1f );

    RealType<ValueType> eps = 0.0001;

    TestDistributions distributions( n );

    DistributionPtr colDist( new NoDistribution( n ) );

    for ( size_t d1 = 0; d1 < distributions.size(); ++d1 )
    {
        DistributionPtr dist1 = distributions[d1];

        if ( dist1->isReplicated() )
        {
            continue;
        }

        for ( size_t d2 = 0; d2 < distributions.size(); ++d2 )
        {

            DistributionPtr dist2 = distributions[d2];

            if ( dist2->isReplicated() )
            {
                continue;
            }

            Matrices<ValueType> allMatrices( context );    // is created by factory

            for ( size_t s = 0; s < allMatrices.size(); ++s )
            {
                Matrix<ValueType>& matrix = *allMatrices[s];

                matrix = cast<ValueType>( csr );

                std::unique_ptr<Matrix<ValueType> > matrix1( matrix.copy() );

                matrix.redistribute( dist1, colDist );

                auto plan = dmemo::redistributePlanByNewDistribution( dist2, dist1 );

                SCAI_LOG_DEBUG( logger, "redistribute plan = " << plan << ", applied to " << matrix )

                matrix.redistribute( plan );
     
                BOOST_CHECK( matrix.isConsistent() );

                matrix1->redistribute( dist2, colDist );

                matrix -= *matrix1;   // sub works local as same distribution

                BOOST_CHECK( matrix.maxNorm() < eps );
            }
        }
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE( disassembleTest )
{
    typedef DefaultReal ValueType;      // assemble/dissamble is not value-type specific

    using namespace hmemo;

    ContextPtr context = Context::getContextPtr();  // test context

    const IndexType m =  8;
    const IndexType n = 15;

    HArray<IndexType> ia(     { 0, 2, 2, 4, 6 } );
    HArray<IndexType> ja(     { 2, 3, 9, 5, 1 } );
    HArray<ValueType> values( { 2, 3, 1, 4, 5 } );

    COOStorage<ValueType> inputCOO( m , n, ia, ja, values );

    TestDistributions rowDists( m );
    TestDistributions colDists( n );

    // Disassemble test for each matrix type with all kind of row and column distributions

    for ( size_t d1 = 0; d1 < rowDists.size(); ++d1 )
    {
        DistributionPtr rowDist = rowDists[d1];

        for ( size_t d2 = 0; d2 < colDists.size(); ++d2 )
        {
            DistributionPtr colDist = colDists[d2];

            Matrices<ValueType> allMatrices( context );    // is created by factory

            for ( size_t s = 0; s < allMatrices.size(); ++s )
            {
                Matrix<ValueType>& matrix = *allMatrices[s];

                matrix.assignDistribute( inputCOO, rowDist, colDist );

                SCAI_LOG_DEBUG( logger, "Disassemble this matrix: " << matrix )

                MatrixAssembly<ValueType> assembly( rowDist->getCommunicatorPtr() );

                matrix.disassemble( assembly );

                COOStorage<ValueType> coo = assembly.buildGlobalCOO( m, n, common::BinaryOp::COPY );

                HArray<IndexType> cooIA = coo.getIA();
                HArray<IndexType> cooJA = coo.getJA();
                HArray<ValueType> cooValues = coo.getValues();

                BOOST_CHECK_EQUAL( cooIA.size(), cooJA.size() );
                BOOST_CHECK_EQUAL( cooIA.size(), cooValues.size() );

                sparsekernel::COOUtils::sort( cooIA, cooJA, cooValues, context );

                BOOST_TEST( hostReadAccess( cooIA ) == hostReadAccess( ia ), per_element() );
                BOOST_TEST( hostReadAccess( cooJA ) == hostReadAccess( ja ), per_element() );
                BOOST_TEST( hostReadAccess( cooValues ) == hostReadAccess( values ), per_element() );
            }
        }
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE_TEMPLATE( hcatTest, ValueType, scai_numeric_test_types )
{
    const IndexType n = 15;
    const IndexType m =  8;

    hmemo::ContextPtr context = hmemo::Context::getContextPtr();  // test context

    common::Math::srandom( 1711 );  // all processors draw same replicated matrix

    auto csr = zero<CSRSparseMatrix<ValueType>>( n, m );

    MatrixCreator::fillRandom( csr, 0.1f );

    CSRSparseMatrix<ValueType> csr2;

    SCAI_LOG_INFO( logger, "hcat( csr, csr ) with csr = " << csr )

    csr2.hcat( csr, csr );

    TestDistributions distributions( n );

    DistributionPtr rowDist( new NoDistribution( n ) );
    DistributionPtr colDist( new NoDistribution( m ) );

    RealType<ValueType> eps = common::TypeTraits<ValueType>::small();

    for ( size_t d1 = 0; d1 < distributions.size(); ++d1 )
    {
        DistributionPtr dist = distributions[d1];

        Matrices<ValueType> allMatrices( context );    // is created by factory

        for ( size_t s = 0; s < allMatrices.size(); ++s )
        {
            Matrix<ValueType>& matrix = *allMatrices[s];

            if ( matrix.getMatrixKind() == MatrixKind::DENSE )
            {
                continue;
            }

            matrix = csr;

            matrix.redistribute( dist, colDist );

            SCAI_LOG_DEBUG( logger, "hcat matrix = " << matrix )

            matrix.hcat( matrix, matrix );

            BOOST_REQUIRE_EQUAL( 2 * n, matrix.getNumRows() );
            BOOST_REQUIRE_EQUAL( m, matrix.getNumColumns() );

            csr2.redistribute( matrix.getRowDistributionPtr(), matrix.getColDistributionPtr() );

            BOOST_CHECK( matrix.maxDiffNorm( csr2 ) < eps );
        }
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE_TEMPLATE( ComplexTest, ValueType, scai_numeric_test_types )
{
    // skip this test if ValueType is not complex as imag would return 0

    if ( !common::isComplex( common::TypeTraits<ValueType>::stype ) )
    {
        return;
    }

    typedef RealType<ValueType> Real;

    dmemo::CommunicatorPtr comm( dmemo::Communicator::getCommunicatorPtr() );
    hmemo::ContextPtr context = hmemo::Context::getContextPtr();  // test context

    const IndexType n = 100;

    Matrices<ValueType> allMatrices( context );    // is created by factory

    dmemo::DistributionPtr dist( new dmemo::BlockDistribution( n, comm ) );

    for ( size_t i = 0; i < allMatrices.size(); ++i )
    {
        Matrix<ValueType>& complexMatrix = *allMatrices[i];

        float fillRate = 0.1f;

        complexMatrix.allocate( dist, dist );
        MatrixCreator::fillRandom( complexMatrix, fillRate  );

        DenseMatrix<Real> x;
        x = real ( complexMatrix );

        DenseMatrix<Real> y;
        y = imag( complexMatrix );

        DenseMatrix<ValueType> z;
        z =  complex( x, y );

        Real diff = complexMatrix.maxDiffNorm( z );

        BOOST_CHECK_EQUAL( diff, 0 );
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE( binaryOpTest )
{
    // only one operation and one ValueType as we here just test correct distributed computing

    typedef DefaultReal ValueType;

    //   Example data
    //
    //    1  0  1   2       2  1  0   2        -1  -1  1 0
    //    2  1  3   0   -   3  0  -1  2    =   -1   1  4  -2
    //    1  0  2   0       1  1  0   0         0   -1  2  0

    const IndexType m = 3;
    const IndexType n = 4;

    HArray<ValueType> data1 ( {  1,  0, 1, 2,  2, 1,  3,  0, 1,  0, 2, 0 } );
    HArray<ValueType> data2 ( {  2,  1, 0, 2,  3, 0, -1,  2, 1,  1, 0, 0 } );
    HArray<ValueType> result( { -1, -1, 1, 0, -1, 1,  4, -2, 0, -1, 2, 0 } );

    DenseStorage<ValueType> storage1( m, n, data1 );
    DenseStorage<ValueType> storage2( m, n, data2 );
    DenseStorage<ValueType> storageR( m, n, result );

    dmemo::CommunicatorPtr comm( dmemo::Communicator::getCommunicatorPtr() );
    hmemo::ContextPtr context = hmemo::Context::getContextPtr();  // test context

    Matrices<ValueType> allMatrices( context );    // is created by factory

    TestDistributions rowDists( m );
    TestDistributions colDists( n );

    auto repRows = std::make_shared<dmemo::NoDistribution>( m );
    auto repCols = std::make_shared<dmemo::NoDistribution>( n );

    for ( size_t i = 0; i < allMatrices.size(); ++i )
    {
        Matrix<ValueType>& matrix1 = *allMatrices[i];

        // for ( size_t d1 = 0; d1 < rowDists.size(); ++d1 )
        for ( size_t d1 = 0; d1 < 1; ++d1 )
        {
            // for ( size_t d2 = 0; d2 < colDists.size(); ++d2 )
            for ( size_t d2 = 0; d2 < 1; ++d2 )
            {
                auto rowDist = rowDists[d1];
                auto colDist = colDists[d2];

                // First test: matrix <binop> SparseMatrix, element-wise

                matrix1.assignDistribute( storage1, rowDist, colDist );

                auto matrix2S = distribute<CSRSparseMatrix<ValueType>>( storage2, rowDist, colDist );

                // Note: matrices for elemen-wise binary op have same distributions

                SCAI_LOG_DEBUG( logger, "binary op, matrix1 = " << matrix1 << ",\nmatrix2 = " << matrix2S )

                matrix1.binaryOp( matrix1, common::BinaryOp::SUB, matrix2S );

                // verify results with replicated data

                matrix1.redistribute( repRows, repCols );

                auto storage = convert<DenseStorage<ValueType>>( matrix1.getLocalStorage() );

                BOOST_TEST( hostReadAccess( storage.getValues() ) == hostReadAccess( storageR.getValues() ), per_element() );

                // second test: matrix <binop> DenseMatrix, element-wise

                matrix1.assignDistribute( storage1, rowDist, colDist );

                auto matrix2D = distribute<DenseMatrix<ValueType>>( storage2, rowDist, colDist );

                // Note: matrices for elemen-wise binary op have same distributions

                matrix1.binaryOp( matrix1, common::BinaryOp::SUB, matrix2D );

                // verify results with replicated data

                matrix1.redistribute( repRows, repCols );

                storage.assign( matrix1.getLocalStorage() );

                BOOST_TEST( hostReadAccess( storage.getValues() ) == hostReadAccess( storageR.getValues() ), per_element() );
            }
        }
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_SUITE_END();
