/**
 * @file NormTest.hpp
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
 * @brief Contains the implementation of the class NormTest
 * @author Alexander Büchel, Micha
 * @date 03.02.2012
 * @since 1.0.0
 */

#include <boost/test/unit_test.hpp>

#include <scai/lama/norm/Norm.hpp>
#include <scai/common/test/TestMacros.hpp>

using namespace scai::lama;

/** Common test class for all derived classes of class Norm.
 *
 *  Each norm has to fulfill some kind of properties, to be correct.
 *  This class provides methods that check for these properties.
 */

static std::string normtestclasses[] =
{ "MaxNormTest", "L1NormTest", "L2NormTest" };

static std::string normtestmethods[] =
{ "positiveHomogeneityTest", "triangleInequalityTest", "ZeroVectorTest" };

class NormTest
{
public:

    NormTest( const Norm& norm )
        : mNorm( norm )
    {
    }

    void positiveHomogeneityTest();
    void triangleInequalityTest();
    void ZeroVectorTest();
    void runTests();

    const Norm& mNorm;
};

#define NORMTEST_COMMONTESTCASES( testinstance )                        \
    {   COMMONTESTCASEINVOKER( testinstance, positiveHomogeneityTest );     \
        COMMONTESTCASEINVOKER( testinstance, triangleInequalityTest );      \
        COMMONTESTCASEINVOKER( testinstance, ZeroVectorTest ); }
