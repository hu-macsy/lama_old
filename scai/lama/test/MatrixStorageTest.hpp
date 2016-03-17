/**
 * @file MatrixStorageTest.hpp
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
 * @brief Contains the implementation of the class MatrixStorageTest
 * @author Thomas Brandes
 * @date 02.03.2012
 * @since 1.0.0
 */

#pragma once

#include <boost/test/unit_test.hpp>

#include <scai/dmemo/Communicator.hpp>
#include <scai/lama/storage/MatrixStorage.hpp>

#include <scai/lama/test/TestMacros.hpp>

static std::string storagetestclasses[] =
{ "CSRStorageTest", "COOStorageTest", "DIAStorageTest", "ELLStorageTest", "JDSStorageTest", "DenseStorageTest" };

static std::string storagetestmethods[] =
{
    "purgeTest", "setCSRDataTest", "buildCSRDataTest", "diagonalTest", "scaleTest", "normTest",
    "vectorMultTest", "jacobiTest", "jacobiHaloTest", "matrixMultTest", "matrixMultTest1", "matrixAddTest"
    "inverseTest", "symmetryTest", "vectorTimesMatrixTest", "numericalTest"
};

/** Test class for MatrixStorage<ValueType>.
 *
 *  LAMA supports different derived classes from MatrixStorage<ValueType> like
 *  CSRStorage<ValueType>, ELLStorage<ValueType>, and so on.
 *
 *  This class tests all common routines for the matrix storages without
 *  taking any knowledge about specific implementation details.
 *
 *  Dependent on the default communicator, it will select between serial
 *  and/or distributed tests.
 *
 * @tparam T is the value type stored in the wrapped container.
 */
template<typename ValueType>
class MatrixStorageTest
{
public:
    /** Constructor of the test.
     *
     *  @param[in] storage is an (derived) object of matrix storage.
     *
     */
    MatrixStorageTest( scai::lama::MatrixStorage<ValueType>& storage )
        : mMatrixStorage( storage )
    {
    }
    ;

    /** Test for MatrixStorage<ValueType>::purge and MatrixStorage<ValueType>::getMemoryUsage */

    void purgeTest();

    /** Test for virtual method MatrixStorage<ValueType>::setCSRData. */

    void setCSRDataTest();

    /** Test for virtual method MatrixStorage<ValueType>::buildCSRData. */

    void buildCSRDataTest();

    /** Test for virtual methods MatrixStorage<ValueType>::setDiagonal and
     *  MatrixStorage<ValueType>::getDiagonal.
     */

    void diagonalTest();

    /** Test for virtual methods MatrixStorage<ValueType>::scale */

    void scaleTest();

    void normTest();

    void vectorMultTest();

    void vectorTimesMatrixTest();

    void numericalTest();

    void jacobiTest();

    void jacobiHaloTest();

    /** simple matrix multiply test for matrices where only diagonals are set. */

    void matrixMultTest();

    /** multiply add test, compare for same results */

    void matrixAddTest();

    /** advanced matrix multiply test, compare for same results by multVector */

    void matrixMultTest1();

    //todo: not implemented --> implement or delete
    //void haloTest( const scai::lama::CommunicatorPtr comm );

    //void replicateTest( const scai::lama::CommunicatorPtr comm );

    //void redistributeTest( const scai::lama::CommunicatorPtr comm );

    void inverseTest();

    void symmetryTest();

    void runTests();

    scai::lama::MatrixStorage<ValueType>& mMatrixStorage;

private:

    /** Test for virtual method MatrixStorage<ValueType>::jacobiIterate.
     *
     *  @param[in] omega is the omega value for the Jacobi iteration.
     *
     *  The omega value allows to test different branches in the
     *  implementation as omega = 0.5 and omega = 1.0 allow for some
     *  kind of optimization.
     */
    void jacobiTest( ValueType omega );

    static void setDenseData( scai::lama::MatrixStorage<ValueType>& storage );

    static void setDenseDataNotSquare( scai::lama::MatrixStorage<ValueType>& storage );

    static void setDenseDataSymmetric( scai::lama::MatrixStorage<ValueType>& storage );

    static void setDenseLocal( scai::lama::MatrixStorage<ValueType>& storage );

    static void setDenseHalo( scai::lama::MatrixStorage<ValueType>& storage );

    static void setDenseRandom( scai::lama::MatrixStorage<ValueType>& storage );

    static void setDenseRandomInverse( scai::lama::MatrixStorage<ValueType>& storage );

    SCAI_LOG_DECL_STATIC_LOGGER( logger );
};

#define MATRIXSTORAGE_COMMONTESTCASES( testinstance )                   \
    {   COMMONTESTCASEINVOKER( testinstance, purgeTest );                   \
        COMMONTESTCASEINVOKER( testinstance, setCSRDataTest );              \
        COMMONTESTCASEINVOKER( testinstance, buildCSRDataTest );            \
        COMMONTESTCASEINVOKER( testinstance, diagonalTest );                \
        COMMONTESTCASEINVOKER( testinstance, scaleTest );                   \
        COMMONTESTCASEINVOKER( testinstance, normTest );                    \
        COMMONTESTCASEINVOKER( testinstance, vectorMultTest );              \
        COMMONTESTCASEINVOKER( testinstance, vectorTimesMatrixTest );       \
        COMMONTESTCASEINVOKER( testinstance, numericalTest );               \
        COMMONTESTCASEINVOKER( testinstance, jacobiTest );                  \
        COMMONTESTCASEINVOKER( testinstance, matrixMultTest );              \
        COMMONTESTCASEINVOKER( testinstance, matrixAddTest );               \
        COMMONTESTCASEINVOKER( testinstance, matrixMultTest1 );             \
        COMMONTESTCASEINVOKER( testinstance, inverseTest ); }
