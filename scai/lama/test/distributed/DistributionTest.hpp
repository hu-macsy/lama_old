/**
 * @file DistributionTest.hpp
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
 * @brief Contains the implementation of the class DistributionTest
 * @author Alexander Büchel
 * @date 30.07.2012
 * @since 1.0.0
 */

#pragma once

#include <boost/test/unit_test.hpp>

#include <scai/dmemo/Distribution.hpp>

static std::string distclasses[] =
{
    "BlockDistributionTest", "GeneralDistributionTest", "NoDistributionTest", "CyclicDistributionTest",
    "GenBlockDistributionTest", "MetisDistributionTest"
};

static std::string distmethods[] =
{ "writeAtTest", "localSizeTest", "global2LocalTest", "printDistributionVector", "local2GlobalTest" };

class DistributionTest
{
public:

    DistributionTest( scai::dmemo::DistributionPtr& dist )
        : mDistributionPtr( dist )
    {
    }

    void writeAtTest();
    void getGlobalSizeTest(); // TODO: missing test ?
    void localSizeTest();
    void local2GlobalTest();
    void global2LocalTest();
    void printDistributionVector();
    void runTests();

private:
    const scai::dmemo::DistributionPtr& mDistributionPtr;
};

#define DISTRIBUTION_COMMONTESTCASES( testinstance )                        \
    {   COMMONTESTCASEINVOKER( testinstance, writeAtTest );                 \
        COMMONTESTCASEINVOKER( testinstance, localSizeTest );               \
        COMMONTESTCASEINVOKER( testinstance, local2GlobalTest );            \
        COMMONTESTCASEINVOKER( testinstance, global2LocalTest );            \
        COMMONTESTCASEINVOKER( testinstance, printDistributionVector ); }
