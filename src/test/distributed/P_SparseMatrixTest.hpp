/**
 * @file P_SparseMatrixTest.hpp
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
 * @brief Contains the implementation of the class P_SparseMatrixTest
 * @author Alexander Büchel
 * @date 10.05.2012
 * @since 1.0.0
 */
#include <boost/test/unit_test.hpp>

#include <lama/matrix/Matrix.hpp>
#include <test/TestMacros.hpp>

using namespace lama;

static std::string psparseMatrixtestclasses[] =
{   "P_CSRSparseMatrixTest", "P_COOSparseMatrixTest", "P_ELLSparseMatrixTest", "P_DIASparseMatrixTest",
    "P_JDSSparseMatrixTest"
};

static std::string psparseMatrixtestmethods[] =
{ "repDistTest", "replicateTest", "assignTest", "transposeTest", "createPoissonTest", "cTorTest" };

template<typename MatrixType>
class P_SparseMatrixTest
{
public:

    P_SparseMatrixTest();

    ~P_SparseMatrixTest();

    void cTorTest();

    void runTests();

    void repDistTest();
    void replicateTest();
    void assignTest();
    void transposeTest();
    void createPoissonTest();

private:
    CommunicatorPtr comm;

protected:
    LAMA_LOG_DECL_STATIC_LOGGER( logger );
};

#define PSPARSEMATRIXTEST_COMMONTESTCASES( testinstance )               \
    {   COMMONTESTCASEINVOKER( testinstance, repDistTest );                 \
        COMMONTESTCASEINVOKER( testinstance, replicateTest );               \
        COMMONTESTCASEINVOKER( testinstance, assignTest );                  \
        COMMONTESTCASEINVOKER( testinstance, transposeTest );               \
        COMMONTESTCASEINVOKER( testinstance, createPoissonTest ); }
