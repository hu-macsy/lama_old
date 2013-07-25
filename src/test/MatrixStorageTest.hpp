/**
 * @file MatrixStorageTest.hpp
 *
 * @license
 * Copyright (c) 2009-2013
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
#include <boost/test/unit_test.hpp>

#include <lama/Communicator.hpp>
#include <lama/storage/MatrixStorage.hpp>

#include <test/TestMacros.hpp>

static std::string storagetestclasses[] =
{ "CSRStorageTest", "COOStorageTest", "DIAStorageTest", "ELLStorageTest", "JDSStorageTest", "DenseStorageTest" };

static std::string storagetestmethods[] =
{   "purgeTest", "emptyTest", "setIdentityTest", "setCSRDataTest", "buildCSRDataTest", "diagonalTest", "scaleTest", "normTest",
    "vectorMultTest", "jacobiTest", "jacobiHaloTest", "matrixMultTest", "matrixMultTest1", "matrixAddTest"
    "writeAtTest", "inverseTest", "vectorTimesMatrixTest"
};

/** Test class for MatrixStorage<T>.
 *
 *  LAMA supports different derived classes from MatrixStorage<T> like
 *  CSRStorage<T>, ELLStorage<T>, and so on.
 *
 *  This class tests all common routines for the matrix storages without
 *  taking any knowledge about specific implementation details.
 *
 *  Dependent on the default communicator, it will select between serial
 *  and/or distributed tests.
 *
 * @tparam T is the value type stored in the wrapped container.
 */
template<typename T>
class MatrixStorageTest
{
public:

    typedef T ValueType; //!< This is the type stored in the wrapped container.

    /** Constructor of the test.
     *
     *  @param[in] storage is an (derived) object of matrix storage.
     *
     */
    MatrixStorageTest( lama::MatrixStorage<T>& storage )
        : mMatrixStorage( storage )
    {
    }
    ;

    /** Test for MatrixStorage<T>::purge and MatrixStorage<T>::getMemoryUsage */

    void purgeTest();

    /** Test for MatrixStorage<T>::hasDiagonalProperty() on empty matrx. */

    void emptyTest();

    /** Test for MatrixStorage<T>::setIdentity. */

    void setIdentityTest();

    /** Test for virtual method MatrixStorage<T>::setCSRData. */

    void setCSRDataTest();

    /** Test for virtual method MatrixStorage<T>::buildCSRData. */

    void buildCSRDataTest();

    /** Test for virtual methods MatrixStorage<T>::setDiagonal and
     *  MatrixStorage<T>::getDiagonal.
     */

    void diagonalTest();

    /** Test for virtual methods MatrixStorage<T>::scale */

    void scaleTest();

    void normTest();

    void vectorMultTest();

    void vectorTimesMatrixTest();

    void jacobiTest();

    void jacobiHaloTest();

    /** simple matrix multiply test for matrices where only diagonals are set. */

    void matrixMultTest();

    /** multiply add test, compare for same results */

    void matrixAddTest();

    /** advanced matrix multiply test, compare for same results by multVector */

    void matrixMultTest1();

    //todo: not implemented --> implement or delete
    //void haloTest( const lama::CommunicatorPtr comm );

    //void replicateTest( const lama::CommunicatorPtr comm );

    //void redistributeTest( const lama::CommunicatorPtr comm );

    void inverseTest();

    void runTests();

    void writeAtTest();

    lama::MatrixStorage<T>& mMatrixStorage;

private:

    /** Test for virtual method MatrixStorage<T>::jacobiIterate.
     *
     *  @param[in] omega is the omega value for the Jacobi iteration.
     *
     *  The omega value allows to test different branches in the
     *  implementation as omega = 0.5 and omega = 1.0 allow for some
     *  kind of optimization.
     */
    void jacobiTest( ValueType omega );

    static void setDenseData( lama::MatrixStorage<T>& storage );

    static void setDenseDataNotSquare( lama::MatrixStorage<T>& storage );

    static void setDenseLocal( lama::MatrixStorage<T>& storage );

    static void setDenseHalo( lama::MatrixStorage<T>& storage );

    static void setDenseRandom( lama::MatrixStorage<T>& storage );

    static void setDenseRandomInverse( lama::MatrixStorage<T>& storage );

    LAMA_LOG_DECL_STATIC_LOGGER( logger );
};

#define MATRIXSTORAGE_COMMONTESTCASES( testinstance )                   \
    {   COMMONTESTCASEINVOKER( testinstance, purgeTest );                   \
        COMMONTESTCASEINVOKER( testinstance, setIdentityTest );             \
        COMMONTESTCASEINVOKER( testinstance, setCSRDataTest );              \
        COMMONTESTCASEINVOKER( testinstance, buildCSRDataTest );            \
        COMMONTESTCASEINVOKER( testinstance, diagonalTest );                \
        COMMONTESTCASEINVOKER( testinstance, scaleTest );                   \
        COMMONTESTCASEINVOKER( testinstance, normTest );                    \
        COMMONTESTCASEINVOKER( testinstance, vectorMultTest );              \
        COMMONTESTCASEINVOKER( testinstance, vectorTimesMatrixTest );       \
        COMMONTESTCASEINVOKER( testinstance, jacobiTest );                  \
        COMMONTESTCASEINVOKER( testinstance, matrixMultTest );              \
        COMMONTESTCASEINVOKER( testinstance, matrixAddTest );               \
        COMMONTESTCASEINVOKER( testinstance, matrixMultTest1 );             \
        COMMONTESTCASEINVOKER( testinstance, writeAtTest );                 \
        COMMONTESTCASEINVOKER( testinstance, inverseTest ); }
