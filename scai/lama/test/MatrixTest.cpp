/**
 * @file MatrixTest.cpp
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
 * @brief Contains the implementation of the class MatrixTest.
 * @author Alexander Büchel
 * @date 01.02.2012
 */

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <scai/common/SCAITypes.hpp>

#include <scai/lama/matrix/ELLSparseMatrix.hpp>
#include <scai/lama/matrix/JDSSparseMatrix.hpp>
#include <scai/lama/matrix/DIASparseMatrix.hpp>
#include <scai/lama/matrix/CSRSparseMatrix.hpp>
#include <scai/lama/matrix/COOSparseMatrix.hpp>
#include <scai/lama/matrix/DenseMatrix.hpp>

#include <scai/lama/expression/MatrixExpressions.hpp>

#include <scai/lama/matutils/MatrixCreator.hpp>

#include <scai/lama/test/TestMacros.hpp>
#include <scai/lama/test/SparseMatrixHelper.hpp>
#include <scai/lama/test/SameMatrixHelper.hpp>
#include <scai/lama/test/TestSparseMatrices.hpp>

#include <cstring>

using namespace scai::lama;
using namespace scai::hmemo;
using namespace scai::dmemo;

BOOST_AUTO_TEST_SUITE( MatrixTest )

/* --------------------------------------------------------------------- */

SCAI_LOG_DEF_LOGGER( logger, "Test.MatrixTest" )

/* --------------------------------------------------------------------- */

typedef SCAI_TEST_TYPE ValueType;

typedef boost::mpl::list<CSRSparseMatrix<ValueType>,
        ELLSparseMatrix<ValueType>,
        COOSparseMatrix<ValueType>,
        JDSSparseMatrix<ValueType>,
        DIASparseMatrix<ValueType>,
        DenseMatrix<ValueType> > MatrixTypes;

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( TypeNameTest, MatrixType, MatrixTypes )
{
    const std::string classTypeName = MatrixType::typeName();
    MatrixType matrix;
    const Matrix* ptrMatrix = &matrix;
    const std::string objTypeName = ptrMatrix->getTypeName();
    SCAI_LOG_INFO( logger, matrix << ": class name = " << classTypeName
                   << ", object name = " << objTypeName );
    BOOST_CHECK( objTypeName == classTypeName );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( RandomTest, MatrixType, MatrixTypes )
{
    typedef typename MatrixType::MatrixValueType ValueType;
    const IndexType numRows = 150;
    const IndexType numColumns = 160;
    const double density = 0.1; // 10% density
    MatrixType matrix( numRows, numColumns );
    MatrixCreator<ValueType>::fillRandom( matrix, density );
    SCAI_LOG_INFO( logger, "Random matrix = " << matrix );
    BOOST_CHECK_EQUAL( numRows, matrix.getNumRows() );
    BOOST_CHECK_EQUAL( numColumns, matrix.getNumColumns() );
    IndexType numValues = matrix.getNumValues();
// should be very unlikely that number of values diverges much from expected number. */
    BOOST_CHECK( density * numRows * numColumns * 2.0 > numValues );
    BOOST_CHECK( density * numRows * numColumns * 0.5 < numValues );
    double sparsityRate = matrix.getSparsityRate();
    BOOST_CHECK( density * 2.0 > sparsityRate );
    BOOST_CHECK( density * 0.5 < sparsityRate );
    size_t memUsage = matrix.getMemoryUsage();
    size_t denseUsage = numRows * numColumns * sizeof( ValueType );

    if ( matrix.getMatrixKind() == Matrix::SPARSE && matrix.getLocalStorage().getFormat() != Format::DIA )
    {
        BOOST_CHECK( memUsage * 2 < denseUsage );
    }
}

/* --------------------------------------------------------------------- */

// Test filename constructor and writing to file for each matrx type
BOOST_AUTO_TEST_CASE_TEMPLATE( ReadWriteTest, MatrixType, MatrixTypes )
{
    CommunicatorPtr comm = Communicator::getCommunicatorPtr(); // get default, MPI or serial
    SCAI_LOG_INFO( logger, "ReadWriteTest for MatrixType = " << typeid( MatrixType ).name() );
    std::string prefix = scai::test::Configuration::getPath();
    SCAI_LOG_INFO( logger, "prefix = " << prefix );
    SCAI_LOG_INFO( logger, "readWriteTest: check loading float matrix" );
// load an available testfile matrix.
    // std::string formattedInputFile = prefix + "/bp__1600.mtx";
    std::string formattedInputFile = prefix + "/impcol_b.mtx";
    MatrixType formattedMatrix( formattedInputFile );
    SCAI_LOG_INFO( logger, "formatted input matrix = " << formattedMatrix );

    /*
     * Its not sure why we need to test with these big matrices here
    if ( formattedMatrix.getMatrixKind() == Matrix::SPARSE )
    {
        // just some additional test, be careful, is too large for dense matrices
        std::string formattedInputFile = prefix + "/2D_poisson_256_formatted.frm";
        std::string xdrInputFile = prefix + "/2D_poisson_256_xdr.frm";
        MatrixType formattedPoisson( formattedInputFile );
        MatrixType xdrPoisson( formattedInputFile );
        testSameMatrix( formattedPoisson, xdrPoisson );
    }
    */

    SCAI_LOG_INFO( logger, "readWriteTest: check writing and loading formatted matrix" );
    std::string formattedFileName = prefix + "/test_matrix_formatted.tmp";

    formattedMatrix.writeToFile( formattedFileName, File::SAMG_FORMAT, scai::common::scalar::FLOAT, scai::common::scalar::INT, scai::common::scalar::INT );
    MatrixType readFormattedMatrix( formattedFileName + ".frm" );
    testSameMatrix( formattedMatrix, readFormattedMatrix );

    SCAI_LOG_INFO( logger, "readWriteTest: check writing and loading binary matrix" );
    std::string binaryFileName = prefix + "/test_matrix_bin.tmp.frm";
// Be careful: binary read must fit to the format that has been used for the write
    formattedMatrix.writeToFile( binaryFileName, File::SAMG_FORMAT, scai::common::scalar::INTERNAL, scai::common::scalar::INT, scai::common::scalar::INT, true );
    MatrixType readBinaryMatrix( binaryFileName );
    testSameMatrix( formattedMatrix, readBinaryMatrix );

    SCAI_LOG_INFO( logger, "readWriteTest: check writing and loading Matrix Market matrix" );
    std::string matrixMarketFileName = prefix + "/test_matrix_mm.tmp";
    formattedMatrix.writeToFile( matrixMarketFileName, File::MATRIX_MARKET, scai::common::scalar::DOUBLE );
    MatrixType readMarketMatrix( matrixMarketFileName + ".mtx" );
    SCAI_LOG_INFO( logger, "mtx matrix: " << readMarketMatrix )
    testSameMatrix( formattedMatrix, readMarketMatrix );

    // remove temporarily created matrix files
    std::remove( (formattedFileName + ".frm").c_str() );
    std::remove( (formattedFileName + ".amg").c_str() );
    std::remove( binaryFileName.c_str() );
    std::remove( (prefix + "/test_matrix_bin.tmp.amg").c_str() );
    std::remove( (matrixMarketFileName + ".mtx" ).c_str() );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( DenseMatrixMultTest, MatrixType, MatrixTypes )
{
    typedef typename MatrixType::MatrixValueType ValueType;
// DenseMatrix = MatrixType * DenseMatrix
    MatrixType matrix1 = TestSparseMatrices::n6m4MatrixE1<ValueType>();
    DenseMatrix<ValueType> matrix2 = TestSparseMatrices::n4m3MatrixE2<ValueType>();
    DenseMatrix<ValueType> matrixP( matrix1 * matrix2 );
    DenseMatrix<ValueType> matrixR = TestSparseMatrices::n6m3MatrixERes<ValueType>();
    testSameMatrix( matrixR, matrixP );
}

/* --------------------------------------------------------------------- */
// TODO Lauri (?) schlaegt fehl.
//BOOST_AUTO_TEST_CASE_TEMPLATE( DenseMatrixAddTest, MatrixType, MatrixTypes )
//{
//
//    typedef typename MatrixType::MatrixValueType ValueType;
//
//    MatrixType matrix1 = TestSparseMatrices::n6m4MatrixE1<ValueType>();
//
//    DenseMatrix<ValueType> matrixP = matrix1 + matrix1;
//
//    DenseMatrix<ValueType> matrixR  = 2.0 * matrix1;
//
//    testSameMatrix( matrixR, matrixP );
//}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END(); /* MatrixTest */
