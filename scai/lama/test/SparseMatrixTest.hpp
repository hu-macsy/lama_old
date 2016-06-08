/**
 * @file SparseMatrixTest.hpp
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
 * @endlicense
 *
 * @brief Contains the implementation of the class SparseMatrixTest
 * @author Alexander Büchel
 * @date 02.03.2012
 */
#include <boost/test/unit_test.hpp>

#include <scai/lama/matrix/Matrix.hpp>
#include <scai/lama/matrix/CSRSparseMatrix.hpp>

#include <scai/lama/test/TestMacros.hpp>

static std::string sparsematrixtestclasses[] =
{
    "CSRSparseMatrixTest", "ELLSparseMatrixTest", "DIASparseMatrixTest", "JDSSParseMatrixTest", "COOSparseMatrixTest",
    "DenseMatrixTest1"
};

static std::string sparsematrixtestmethods[] =
{
    "clearTest", "cTorTest", "testConversions", "testMultiplication", "MatrixExpressionTest", "MatrixCtorExpressionTest",
    "writeAtTest", "scaleTest"
};

/** Test routines for SparseMatrix.
 *
 *  All tests run also fine for DenseMatrix.
 *
 *  @todo: rename SparseMatrixTest to MatrixTest
 *  @todo: remove routines of DenseMatrixTest tested here
 *  @todo: DenseMatrixTest should work also with SparseMatrix (at least most tests)
 *  @todo: identify tests either only for SparseMatrix or only for DenseMatrix
 *  @todo: test also with CUDAContext
 *  @todo: test also with MPICommunicator
 *
 *  @tparam MatrixType any full class derived from Matrix
 *
 */

template<typename MatrixType>
class SparseMatrixTest
{
public:
    SparseMatrixTest( MatrixType& matrix )
        : mMatrix( matrix )
    {
    }

    //void readWriteTest(); -> not implemented?
    void testMultiplication();
    void clearTest();
    void cTorTest();
    void testConversions();
    void MatrixExpressionTest();
    void MatrixCtorExpressionTest();
    //todo: not implemented --> implement or delete
    //void testSortRows();
    void writeAtTest();
    void scaleTest();

    void runTests();

    void setUp();

    template<typename mt> void testConversionsImpl();
    void matrixMultTestImpl( const scai::lama::Matrix& a, const scai::lama::Matrix& b, const scai::lama::Matrix& result );
    void matrixEqualityCheck( const MatrixType& a, const MatrixType& b );

    //bool testDiagonalProperty( const MatrixType& a);

    MatrixType& mMatrix;
    std::string m_FormattedInputFile;
    std::string m_XDRInputFile;
    std::string m_TestOutputFileFor;
    std::string m_TestOutputFileUnf;
    std::string m_TestOutputFileXDR;
    std::string m_TestOutputFileMM;

protected:

    SCAI_LOG_DECL_STATIC_LOGGER( logger )
};

#define SPARSEMATRIX_COMMONTESTCASES( testinstance )                        \
    {   COMMONTESTCASEINVOKER( testinstance, clearTest );                   \
        COMMONTESTCASEINVOKER( testinstance, cTorTest );                    \
        COMMONTESTCASEINVOKER( testinstance, testMultiplication );          \
        COMMONTESTCASEINVOKER( testinstance, MatrixExpressionTest );        \
        COMMONTESTCASEINVOKER( testinstance, MatrixCtorExpressionTest );    \
        COMMONTESTCASEINVOKER( testinstance, writeAtTest );                 \
        COMMONTESTCASEINVOKER( testinstance, scaleTest );                   \
        COMMONTESTCASEINVOKER( testinstance, testConversions ); }
