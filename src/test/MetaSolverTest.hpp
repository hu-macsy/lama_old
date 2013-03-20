/**
 * @file MetaSolverTest.hpp
 *
 * @license
 * Copyright (c) 2012
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
 * @brief Contains the implementation of the class JacobiTest.cpp
 * @author: Kai Buschulte
 * @date 07.05.2012
 * $
 **/

#include <boost/test/unit_test.hpp>

#include <lama/solver/MetaSolver.hpp>

#include <test/TestMacros.hpp>

namespace MetaSolverTest
{

using namespace boost;
using namespace lama;

/* --------------------------------------------------------------------- */

class MetaSolverTestImpl
{
public:
    // double only, because float is not supported by SimpleAMG.
    // Different ValueTypes should be covered by the native Tests.
    typedef double ValueType;

    MetaSolverTestImpl( std::string& config, SolverPtr native );
    virtual ~MetaSolverTestImpl();

    void configTest();
    void initializeTest();
    void solveTest();
    void testCompareNative();
    void writeAtTest();

    void runTests();

protected:
    std::string mConfiguration;
    SolverPtr mNative;
    MatrixPtr mMatrix;
    VectorPtr mSolution;
    VectorPtr mRhs;
    SolverPtr mMetaSolver;
};

} // namespace MetaSolverTestSuite
