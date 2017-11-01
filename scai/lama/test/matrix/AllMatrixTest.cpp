/**
 * @file AllMatrixTest.cpp
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
 * @brief Test cases applied to each matrix class, i.e. test (virtual) methods of _Matrix
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
#include <scai/dmemo/Redistributor.hpp>
#include <scai/dmemo/test/TestDistributions.hpp>

#include <scai/common/TypeTraits.hpp>
#include <scai/common/Math.hpp>

#include <scai/lama/storage/DenseStorage.hpp>
#include <scai/lama/Scalar.hpp>
#include <scai/lama/expression/all.hpp>
#include <scai/lama/matrix/CSRSparseMatrix.hpp>
#include <scai/lama/matrix/DIASparseMatrix.hpp>
#include <scai/lama/matutils/MatrixCreator.hpp>

#include <scai/logging.hpp>


using namespace scai;
using namespace lama;
using namespace dmemo;
using utilskernel::LArray;

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
    _Matrices allMatrices;    // is created by factory
    size_t nFormats = static_cast<size_t>( Format::UNDEFINED );
    size_t nTypes   = SCAI_COMMON_COUNT_NARG( SCAI_NUMERIC_TYPES_HOST );
    nFormats--;   // SPARSE_ASSEMBLY_STORAGE not used for a matrix
    nFormats--;   // STENCIL_STORAGE not used for a matrix
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

BOOST_AUTO_TEST_CASE( _copyTest )
{
    hmemo::ContextPtr context = hmemo::Context::getContextPtr();  // test context
    // For copy we just take one arithmetic type to reduce number of test cases
    common::scalar::ScalarType stype = common::TypeTraits<SCAI_TEST_TYPE>::stype;
    _Matrices allMatrices( stype, context );    // is created by factory
    SCAI_LOG_INFO( logger, "Test " << allMatrices.size() << "  matrices for copy method" )

    for ( size_t s = 0; s < allMatrices.size(); ++s )
    {
        _Matrix& matrix = *allMatrices[s];
        initMatrix( matrix, "BLOCK", "BLOCK" );
        MatrixPtr copyMatrix( matrix.copy() );
        SCAI_LOG_DEBUG( logger, "copyTest: " << matrix << " with copy " << *copyMatrix );
        // verify for same matrix
        BOOST_CHECK_EQUAL( matrix.getRowDistribution(), copyMatrix->getRowDistribution() );
        BOOST_CHECK_EQUAL( matrix.getColDistribution(), copyMatrix->getColDistribution() );
        BOOST_CHECK_EQUAL( matrix.getFormat(), copyMatrix->getFormat() );
        BOOST_CHECK_EQUAL( matrix.getValueType(), copyMatrix->getValueType() );
        BOOST_CHECK_EQUAL( matrix.getContextPtr(), copyMatrix->getContextPtr() );
        Scalar maxDiff = copyMatrix->_maxDiffNorm( matrix );
        // value should be exactly 0
        BOOST_CHECK_EQUAL( 0, maxDiff.getValue<float>() );
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
        typename Matrix<ValueType>::Ptr copyMatrix( matrix.copy() );
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

BOOST_AUTO_TEST_CASE( l1NormTest )
{
    const IndexType N = 8;
    const Scalar scale( 2 );
    const Scalar expectedNorm = scale * Scalar( N );
    hmemo::ContextPtr context = hmemo::Context::getContextPtr();  // test context
    _Matrices allMatrices( context );    // is created by factory
    SCAI_LOG_INFO( logger, "Test " << allMatrices.size() << "  matrices for l1Norm" )

    for ( size_t s = 0; s < allMatrices.size(); ++s )
    {
        _Matrix& matrix = *allMatrices[s];
        matrix.setIdentity( N );
        matrix *= scale;
        SCAI_LOG_DEBUG( logger, "Test l1Norm for this matrix: " << matrix )
        Scalar l1Norm = matrix._l1Norm();
        SCAI_CHECK_CLOSE( expectedNorm, l1Norm, 0.001 );
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE_TEMPLATE( l2NormTest, ValueType, scai_numeric_test_types )
{
    typedef typename Matrix<ValueType>::RealType RealType;

    const IndexType N = 8;
    const ValueType scale( 2 );
    const RealType expectedNorm = common::Math::sqrt( N * scale * scale );

    hmemo::ContextPtr context = hmemo::Context::getContextPtr();  // test context

    Matrices<ValueType> allMatrices( context );    // is created by factory

    SCAI_LOG_DEBUG( logger, "Test " << allMatrices.size() << "  matrices for l2Norm" )

    for ( size_t s = 0; s < allMatrices.size(); ++s )
    {
        Matrix<ValueType>& matrix = *allMatrices[s];
        matrix.setIdentity( N );
        matrix *= scale;
        SCAI_LOG_DEBUG( logger, "Test l2Norm for this matrix: " << matrix )
        RealType l2Norm = matrix.l2Norm();
        BOOST_CHECK_CLOSE( expectedNorm, l2Norm, 0.001 );
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE( maxNormTest )
{
    const IndexType N = 8;
    const Scalar scale( 2 );
    const Scalar expectedNorm = scale;
    hmemo::ContextPtr context = hmemo::Context::getContextPtr();  // test context
    _Matrices allMatrices( context );    // is created by factory
    SCAI_LOG_INFO( logger, "Test " << allMatrices.size() << "  matrices for maxNorm" )

    for ( size_t s = 0; s < allMatrices.size(); ++s )
    {
        _Matrix& matrix = *allMatrices[s];
        matrix.setIdentity( N );
        matrix *= scale;
        SCAI_LOG_DEBUG( logger, "Test maxNorm for this matrix: " << matrix )
        Scalar maxNorm = matrix._maxNorm();
        SCAI_CHECK_CLOSE( expectedNorm, maxNorm, 0.001 );
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

            matrixT.assignTranspose( matrix );
            SCAI_LOG_DEBUG( logger, "transposeTest: " << matrixT << " , orig is " << matrix );
            BOOST_CHECK_EQUAL( matrix.getRowDistribution(), matrixT.getColDistribution() );
            BOOST_CHECK_EQUAL( matrix.getColDistribution(), matrixT.getRowDistribution() );
            Matrix<ValueType>::Ptr matrixTT( matrix.newMatrix() );
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
        Matrix<ValueType>::Ptr copyMatrix( matrix.copy() );
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

        typename Matrix<ValueType>::Ptr matrix2Ptr( matrix1.copy() );
        typename Matrix<ValueType>::Ptr matrix3Ptr( matrix1.newMatrix() );

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

BOOST_AUTO_TEST_CASE( assignMultTest )
{
    hmemo::ContextPtr context = hmemo::Context::getContextPtr();  // test context

    common::scalar::ScalarType stype = common::TypeTraits<SCAI_TEST_TYPE>::stype;

    _Matrices allMatrices( stype, context );    // is created by factory

    SCAI_LOG_INFO( logger, "Test " << allMatrices.size() << "  matrices for assign operator tests" )

    for ( size_t s = 0; s < allMatrices.size(); ++s )
    {
        _Matrix& matrix1 = *allMatrices[s];

        if ( matrix1.getMatrixKind() == MatrixKind::DENSE )
        {
            continue;
        }

        initMatrix( matrix1, "BLOCK", "NO" );

        CSRSparseMatrix<SCAI_TEST_TYPE> unityLeft;
        CSRSparseMatrix<SCAI_TEST_TYPE> unityRight;

        unityLeft.setIdentity( matrix1.getRowDistributionPtr() );
        unityRight.setIdentity( matrix1.getColDistributionPtr() );

        MatrixPtr matrix2Ptr( matrix1.newMatrix() );

        _Matrix& matrix2 = *matrix2Ptr;

        matrix2 = matrix1 * unityRight;   // not for Dense

        BOOST_CHECK_EQUAL( matrix2.getRowDistributionPtr(), matrix1.getRowDistributionPtr() );

        matrix2 = matrix1 * unityRight + matrix1;

        BOOST_CHECK_EQUAL( matrix2.getRowDistributionPtr(), matrix1.getRowDistributionPtr() );
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE( checkSymmetryTest )
{
    hmemo::ContextPtr context = hmemo::Context::getContextPtr();  // test context

    _Matrices allMatrices( context );    // is created by factory

    SCAI_LOG_INFO( logger, "Test " << allMatrices.size() << "  matrices for checkSymmetry" )

    for ( size_t s = 0; s < allMatrices.size(); ++s )
    {
        _Matrix& matrix = *allMatrices[s];

        matrix.setIdentity( 5 );
        BOOST_CHECK( matrix.checkSymmetry() );

        MatrixCreator::buildPoisson2D( matrix, 5, 3, 3 );
        BOOST_CHECK( matrix.checkSymmetry() );
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE( setDiagonalPropertyTest )
{
    const IndexType n1 = 3;
    const IndexType n2 = 4;

    hmemo::ContextPtr context = hmemo::Context::getContextPtr();  // test context

    _Matrices allMatrices( context );    // is created by factory

    TestDistributions testDistributions( n1 * n2 );
    DistributionPtr repDist( new NoDistribution( n1 * n2 ) );

    SCAI_LOG_INFO( logger, "Test " << allMatrices.size() << "  matrices for checkSymmetry" )

    for ( size_t s = 0; s < allMatrices.size(); ++s )
    {
        _Matrix& matrix = *allMatrices[s];

        if ( matrix.getMatrixKind() == MatrixKind::DENSE )
        {
            continue;   // Dense does not support first column indexes
        }

        if ( matrix.getFormat() == Format::DIA )
        {
            continue;   // DIA does not support first column indexes
        }

        for ( size_t i = 0; i < testDistributions.size(); ++i )
        {
            DistributionPtr dist = testDistributions[i];

            matrix.clear();

            MatrixCreator::buildPoisson2D( matrix, 5, n1, n2 );

            matrix.setDiagonalProperty();

            matrix.redistribute( dist, repDist );

            utilskernel::LArray<IndexType> myGlobalIndexes1;
            utilskernel::LArray<IndexType> myGlobalIndexes2;

            matrix.getLocalStorage().getFirstColumnIndexes( myGlobalIndexes1 );
            dist->getOwnedIndexes( myGlobalIndexes2 );

            BOOST_CHECK_EQUAL( myGlobalIndexes1.size(), myGlobalIndexes2.size() );
            BOOST_CHECK_EQUAL( IndexType( 0 ), myGlobalIndexes1.maxDiffNorm( myGlobalIndexes2 ) );
        }
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE( diagonalTest )
{
    const IndexType n1 = 3;
    const IndexType n2 = 4;

    hmemo::ContextPtr context = hmemo::Context::getContextPtr();  // test context

    _Matrices allMatrices( context );    // is created by factory

    TestDistributions testDistributions( n1 * n2 );

    SCAI_LOG_INFO( logger, "Test " << allMatrices.size() << "  matrices for checkSymmetry" )

    for ( size_t s = 0; s < allMatrices.size(); ++s )
    {
        _Matrix& matrix = *allMatrices[s];

        if ( matrix.getFormat() == Format::DIA )
        {
            continue;  // DIA has problems with diagonal property
        }

        matrix.setCommunicationKind( SyncKind::SYNCHRONOUS );

        for ( size_t i = 0; i < testDistributions.size(); ++i )
        {
            DistributionPtr dist = testDistributions[i];

            matrix.clear();

            MatrixCreator::buildPoisson2D( matrix, 5, n1, n2 );

            matrix.redistribute( dist, dist );

            SCAI_LOG_DEBUG( logger, "diagonalTest for " << matrix )

            VectorPtr xPtr ( _Vector::getVector( VectorKind::DENSE, matrix.getValueType() ) );
            VectorPtr y1Ptr( _Vector::getVector( VectorKind::DENSE, matrix.getValueType() ) );
            VectorPtr y2Ptr( _Vector::getVector( VectorKind::DENSE, matrix.getValueType() ) );
            VectorPtr dPtr ( _Vector::getVector( VectorKind::DENSE, matrix.getValueType() ) );

            _Vector& x  = *xPtr;
            _Vector& y1 = *y1Ptr;
            _Vector& y2 = *y2Ptr;
            _Vector& d  = *dPtr;

            x.allocate( matrix.getColDistributionPtr() );
            x  = Scalar( 1 );

            y1 = matrix * x;

            matrix.getDiagonal( d );

            BOOST_CHECK_EQUAL( d.getDistribution(), matrix.getRowDistribution() );

            matrix.setDiagonal( 0 );

            // Now we can prove y2 = matrix * x + diagonal must be same as y1

            y2 = matrix * x + d;
            y2 = y2 - y1;
            BOOST_CHECK( y2._maxNorm() < Scalar( 0.0001 ) );

            // Write back modified diagonal and check result

            d += Scalar( 1 );

            matrix.setDiagonal( d );

            y2 = matrix * x - y1;
            y2 += -1;

            BOOST_CHECK( y2._maxNorm() < Scalar( 0.0001 ) );
        }
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE( getRowTest )
{
    const IndexType nRows = 10;
    const IndexType nCols = 8;

    hmemo::ContextPtr context = hmemo::Context::getContextPtr();  // test context

    CSRSparseMatrix<RealType> csr( nRows, nCols );
    MatrixCreator::fillRandom( csr, 0.1f );

    TestDistributions rowDistributions( nRows );
    TestDistributions colDistributions( nCols );

    for ( size_t i = 0; i < rowDistributions.size(); ++i )
    {
        DistributionPtr rowDist = rowDistributions[i];

        for ( size_t j = 0; j < colDistributions.size(); ++j )
        {
            DistributionPtr colDist = colDistributions[j];

            _Matrices allMatrices( context );    // is created by factory

            for ( size_t s = 0; s < allMatrices.size(); ++s )
            {
                _Matrix& matrix = *allMatrices[s];

                matrix = csr;

                matrix.redistribute( rowDist, colDist );

                SCAI_LOG_INFO( logger, "getRowTest for this matrix: " << matrix )

                // get each row and subtract it

                VectorPtr row ( _Vector::getVector( VectorKind::DENSE, matrix.getValueType() ) );

                for ( IndexType iRow = 0; iRow < matrix.getNumRows(); ++iRow )
                {
                    matrix.getRow( *row, iRow );
                    // BOOST_REQUIRE_EQUAL( nCols, row->size() );
                    // BOOST_REQUIRE( row->isConsistent() );
                    matrix.setRow( *row, iRow, common::binary::SUB );
                }

                // the final matrix should be zero

                BOOST_CHECK( matrix._maxNorm() < Scalar( 0.001 ) );
            }
        }
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE( reduceTest )
{
    const IndexType nRows = 10;
    const IndexType nCols = 8;

    hmemo::ContextPtr context = hmemo::Context::getContextPtr();  // test context

    CSRSparseMatrix<RealType> csr( nRows, nCols );
    MatrixCreator::fillRandom( csr, 0.1f );

    common::binary::BinaryOp reduceOp = common::binary::ADD;
    common::unary::UnaryOp   elemOp   = common::unary::SQR;

    TestDistributions rowDistributions( nRows );
    TestDistributions colDistributions( nCols );

    for ( IndexType dim = 0; dim < 2; ++dim )
    {
        DenseVector<RealType> sRow;
        csr.reduce( sRow, dim, reduceOp, elemOp );

        for ( size_t i = 0; i < rowDistributions.size(); ++i )
        {
            DistributionPtr rowDist = rowDistributions[i];

            for ( size_t j = 0; j < colDistributions.size(); ++j )
            {
                DistributionPtr colDist = colDistributions[j];

                _Matrices allMatrices( context );    // is created by factory

                for ( size_t s = 0; s < allMatrices.size(); ++s )
                {
                    _Matrix& matrix = *allMatrices[s];

                    matrix = csr;

                    matrix.redistribute( rowDist, colDist );

                    VectorPtr row ( _Vector::getVector( VectorKind::DENSE, matrix.getValueType() ) );

                    // reduce on the parallel matrix

                    matrix.reduce( *row, dim, reduceOp, elemOp );

                    if ( dim == 0 )
                    {
                        BOOST_CHECK_EQUAL( row->size(), nRows );
                    }
                    else
                    {
                        BOOST_CHECK_EQUAL( row->size(), nCols );
                    }

                    row->redistribute( sRow.getDistributionPtr() );

                    if ( ! ( row->_maxDiffNorm( sRow ) < Scalar( 0.001 ) ) )
                    {
                        SCAI_LOG_ERROR( logger,  "i = " << i << ", j = " << j << ", s = " << s
                                        << ", fail for this matrix: " << matrix << std::endl );
                    }

                    BOOST_CHECK( row->_maxDiffNorm( sRow ) < Scalar( 0.001 ) );
                }
            }
        }
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE( getColTest )
{
    const IndexType nRows = 10;
    const IndexType nCols = 8;

    hmemo::ContextPtr context = hmemo::Context::getContextPtr();  // test context

    CSRSparseMatrix<RealType> csr( nRows, nCols );
    MatrixCreator::fillRandom( csr, 0.1f );

    TestDistributions rowDistributions( nRows );
    TestDistributions colDistributions( nCols );

    for ( size_t i = 0; i < rowDistributions.size(); ++i )
    {
        DistributionPtr rowDist = rowDistributions[i];

        for ( size_t j = 0; j < colDistributions.size(); ++j )
        {
            DistributionPtr colDist = colDistributions[j];

            _Matrices allMatrices( context );    // is created by factory

            for ( size_t s = 0; s < allMatrices.size(); ++s )
            {
                _Matrix& matrix = *allMatrices[s];

                if ( matrix.getMatrixKind() != MatrixKind::DENSE )
                {
                    // continue;
                }

                matrix = csr;

                matrix.redistribute( rowDist, colDist );

                // get each row and subtract it

                VectorPtr col ( _Vector::getVector( VectorKind::DENSE, matrix.getValueType() ) );

                SCAI_LOG_INFO( logger, "getColTest for this matrix: " << matrix )

                for ( IndexType iCol = 0; iCol < matrix.getNumColumns(); ++iCol )
                {
                    matrix.getColumn( *col, iCol );
                    matrix.setColumn( *col, iCol, common::binary::SUB );
                }

                // the final matrix should be zero

                BOOST_CHECK( matrix._maxNorm() < Scalar( 0.001 ) );
            }
        }
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE( getTest )
{
    const IndexType nRows = 10;
    const IndexType nCols = 8;

    hmemo::ContextPtr context = hmemo::Context::getContextPtr();  // test context

    CSRSparseMatrix<RealType> csr( nRows, nCols );
    MatrixCreator::fillRandom( csr, 0.1f );

    TestDistributions rowDistributions( nRows );
    TestDistributions colDistributions( nCols );

    for ( size_t rd = 0; rd < rowDistributions.size(); ++rd )
    {
        DistributionPtr rowDist = rowDistributions[rd];

        for ( size_t cd = 0; cd < colDistributions.size(); ++cd )
        {
            DistributionPtr colDist = colDistributions[cd];

            _Matrices allMatrices( context );    // is created by factory

            for ( size_t s = 0; s < allMatrices.size(); ++s )
            {
                _Matrix& matrix = *allMatrices[s];

                matrix = csr;

                matrix.redistribute( rowDist, colDist );

                SCAI_LOG_DEBUG( logger, "getTest for this matrix: " << matrix << ", max = " << matrix._maxNorm() )

                // get each element and subract it

                for ( IndexType iRow = 0; iRow < matrix.getNumRows(); ++iRow )
                {
                    for ( IndexType jCol = 0; jCol < matrix.getNumColumns(); ++jCol )
                    {
                        Scalar s1 = matrix.getValue( iRow, jCol );
                        Scalar s2 = csr.getValue( iRow, jCol );

                        Scalar diff = abs( s1 - s2 );

                        BOOST_CHECK( diff < Scalar( 0.001 ) );
                    }
                }
            }
        }
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE( getSetTest )
{
    const IndexType nRows = 10;
    const IndexType nCols = 8;

    hmemo::ContextPtr context = hmemo::Context::getContextPtr();  // test context

    CSRSparseMatrix<RealType> csr( nRows, nCols );
    MatrixCreator::fillRandom( csr, 0.1f );

    TestDistributions rowDistributions( nRows );
    TestDistributions colDistributions( nCols );

    for ( size_t rd = 0; rd < rowDistributions.size(); ++rd )
    {
        DistributionPtr rowDist = rowDistributions[rd];

        for ( size_t cd = 0; cd < colDistributions.size(); ++cd )
        {
            DistributionPtr colDist = colDistributions[cd];

            _Matrices allMatrices( context );    // is created by factory

            for ( size_t s = 0; s < allMatrices.size(); ++s )
            {
                _Matrix& matrix = *allMatrices[s];

                matrix = csr;

                if ( matrix.getMatrixKind() == MatrixKind::SPARSE )
                {
                    continue;
                }

                if ( matrix.getValueType() != common::scalar::FLOAT )
                {
                    continue;
                }

                matrix.redistribute( rowDist, colDist );

                SCAI_LOG_DEBUG( logger, "getSetTest for this matrix: " << matrix << ", max = " << matrix._maxNorm() )

                // get each element and subract it

                for ( IndexType iRow = 0; iRow < matrix.getNumRows(); ++iRow )
                {
                    for ( IndexType jCol = 0; jCol < matrix.getNumColumns(); ++jCol )
                    {
                        Scalar s = matrix.getValue( iRow, jCol );

                        if ( s != Scalar( 0 ) )
                        {
                            matrix.setValue( iRow, jCol, s, common::binary::SUB );
                        }
                    }
                }

                // the final matrix should be zero

                BOOST_CHECK( matrix._maxNorm() < Scalar( 0.001 ) );
            }
        }
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE( setCSRDataTest )
{
    hmemo::ContextPtr context = hmemo::Context::getContextPtr();  // test context

    typedef RealType ValueType;

    const IndexType nRows = 15;
    const IndexType nCols = 8;

    std::srand( 1311 );

    CSRSparseMatrix<ValueType> csr( nRows, nCols );
    MatrixCreator::fillRandom( csr, 0.1f );

    dmemo::DistributionPtr colDist( new NoDistribution( nCols ) );

    TestDistributions testDistributions( nRows );

    for ( size_t i = 0; i < testDistributions.size(); ++i )
    {
        DistributionPtr dist = testDistributions[i];

        csr.redistribute( dist, csr.getColDistributionPtr() );

        CSRStorage<ValueType> localCSR = csr.getLocalStorage();

        const hmemo::HArray<IndexType>& ia = localCSR.getIA();
        const hmemo::HArray<IndexType>& ja = localCSR.getJA();
        const hmemo::HArray<ValueType>& values = localCSR.getValues();
        const IndexType numValues = localCSR.getNumValues();

        // now we have local CSR data to set

        _Matrices allMatrices( context );    // is created by factory

        for ( size_t k = 0; k < allMatrices.size(); ++k )
        {
            _Matrix& mat = *allMatrices[ k ];

            mat.setCSRData( dist, colDist, numValues, ia, ja, values );

            BOOST_CHECK( mat.isConsistent() );
        }
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE( setDIADataTest )
{
    hmemo::ContextPtr context = hmemo::Context::getContextPtr();  // test context

    typedef RealType ValueType;

    const IndexType nRows = 15;
    const IndexType nCols = 8;

    std::srand( 1311 );

    DIASparseMatrix<ValueType> dia( nRows, nCols );
    MatrixCreator::fillRandom( dia, 0.1f );

    dmemo::DistributionPtr colDist( new NoDistribution( nCols ) );

    TestDistributions testDistributions( nRows );

    for ( size_t i = 0; i < testDistributions.size(); ++i )
    {
        DistributionPtr dist = testDistributions[i];

        dia.redistribute( dist, dia.getColDistributionPtr() );

        DIAStorage<ValueType> localDIA = dia.getLocalStorage();

        SCAI_LOG_INFO( logger, "Local DIAData: " << localDIA )

        const IndexType numDiagonals = localDIA.getNumDiagonals();
        const hmemo::HArray<IndexType>& offsets = localDIA.getOffsets();
        const hmemo::HArray<ValueType>& values = localDIA.getValues();

        // now we have local DIA data to set

        _Matrices allMatrices( context );    // is created by factory

        for ( size_t k = 0; k < allMatrices.size(); ++k )
        {
            _Matrix& mat = *allMatrices[ k ];

            if ( mat.getMatrixKind() != MatrixKind::DENSE )
            {
                continue;
            }

            SCAI_LOG_INFO( logger, "setDIAData for this mat: " << mat )

            mat.setDIAData( dist, colDist, numDiagonals, offsets, values );

            BOOST_CHECK( mat.isConsistent() );
        }
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE( redistributeTest )
{
    const IndexType n = 15;   // only square matrices here

    hmemo::ContextPtr context = hmemo::Context::getContextPtr();  // test context

    CSRSparseMatrix<RealType> csr( n, n );
    MatrixCreator::fillRandom( csr, 0.1f );

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

            _Matrices allMatrices( context );    // is created by factory

            for ( size_t s = 0; s < allMatrices.size(); ++s )
            {
                _Matrix& matrix = *allMatrices[s];

                matrix = csr;

                common::unique_ptr<_Matrix> matrix1( matrix.copy() );

                matrix.redistribute( dist1, colDist);

                dmemo::Redistributor redistributor( dist2, dist1 );

                SCAI_LOG_DEBUG( logger, "redistributor = " << redistributor << ", applied to " << matrix )

                matrix.redistribute( redistributor );
 
                matrix1->redistribute( dist2, colDist );

                matrix -= *matrix1;   // sub works local as same distribution

                BOOST_CHECK( matrix._maxNorm() < 0.001 );
            }
        }
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE( hcatTest )
{
    return;

    const IndexType n = 15; 
    const IndexType m =  8;  

    hmemo::ContextPtr context = hmemo::Context::getContextPtr();  // test context

    common::Math::srandom( 1711 );  // all processors draw same replicated matrix

    CSRSparseMatrix<RealType> csr( n, m );
    MatrixCreator::fillRandom( csr, 0.1f );

    CSRSparseMatrix<RealType> csr2;
    csr2.hcat( csr, csr );           

    TestDistributions distributions( n );

    DistributionPtr rowDist( new NoDistribution( n ) );
    DistributionPtr colDist( new NoDistribution( m ) );

    for ( size_t d1 = 0; d1 < distributions.size(); ++d1 )
    {
        DistributionPtr dist = distributions[d1];

        _Matrices allMatrices( context );    // is created by factory

        for ( size_t s = 0; s < allMatrices.size(); ++s )
        {
            _Matrix& matrix = *allMatrices[s];

            if ( matrix.getMatrixKind() == MatrixKind::DENSE )
            {
                continue;
            }

            matrix = csr;

            matrix.redistribute( dist, colDist);

            matrix.hcat( matrix, matrix );

            BOOST_REQUIRE_EQUAL( 2 * n, matrix.getNumRows() );
            BOOST_REQUIRE_EQUAL( m, matrix.getNumColumns() );

            csr2.redistribute( matrix.getRowDistributionPtr(), matrix.getColDistributionPtr() );

            BOOST_CHECK( matrix._maxDiffNorm( csr2 ) < 0.001 );
        }
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_SUITE_END();
