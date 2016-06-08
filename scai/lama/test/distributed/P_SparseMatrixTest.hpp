/**
 * @file P_SparseMatrixTest.hpp
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
 * @brief Contains the implementation of the class P_SparseMatrixTest
 * @author Alexander Büchel
 * @date 10.05.2012
 */
#include <boost/test/unit_test.hpp>

#include <scai/lama/matrix/Matrix.hpp>
#include <scai/lama/test/TestMacros.hpp>

static std::string psparseMatrixtestclasses[] =
{
    "P_CSRSparseMatrixTest", "P_COOSparseMatrixTest", "P_ELLSparseMatrixTest", "P_DIASparseMatrixTest",
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
    scai::dmemo::CommunicatorPtr comm;

protected:
    SCAI_LOG_DECL_STATIC_LOGGER( logger );
};

#define PSPARSEMATRIXTEST_COMMONTESTCASES( testinstance )                   \
    {   COMMONTESTCASEINVOKER( testinstance, repDistTest );                 \
        COMMONTESTCASEINVOKER( testinstance, replicateTest );               \
        COMMONTESTCASEINVOKER( testinstance, assignTest );                  \
        COMMONTESTCASEINVOKER( testinstance, transposeTest );               \
        COMMONTESTCASEINVOKER( testinstance, createPoissonTest ); }
