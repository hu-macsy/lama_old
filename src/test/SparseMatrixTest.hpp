/**
 * @file SparseMatrixTest.hpp
 *
 * @license
 * Copyright (c) 2009-2015
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 * @endlicense
 *
 * @brief Contains the implementation of the class SparseMatrixTest
 * @author Alexander Büchel
 * @date 02.03.2012
 * @since 1.0.0
 */
#include <boost/test/unit_test.hpp>

#include <lama/matrix/Matrix.hpp>
#include <lama/matrix/CSRSparseMatrix.hpp>

#include <test/TestMacros.hpp>

using namespace lama;

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
    void matrixMultTestImpl( const Matrix& a, const Matrix& b, const Matrix& result );
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

    LAMA_LOG_DECL_STATIC_LOGGER( logger )
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
