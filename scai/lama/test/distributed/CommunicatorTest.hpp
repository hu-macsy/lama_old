/**
 * @file CommunicatorTest.hpp
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
 * @brief Contains the implementation of the class CommunicatorTest
 * @author Alexander Büchel
 * @date 09.05.2012
 * @since 1.0.0
 */

#include <boost/test/unit_test.hpp>

#include <scai/dmemo/Communicator.hpp>

/* --------------------------------------------------------------------- */

static std::string commtestclasses[] = { "NoCommunicatorTest" };
static std::string commtestmethods[] =
{
    "swapTest", "gatherTest", "gatherVTest", "scatterTest", "scatterVTest", "bcastTest", "shiftASyncTest", "shiftTest",
    "updateHaloTest", "bcastStringTest", "buildHaloTest", "allocatePlanTest", "computeOwnersTest", "CommunicatorCtorTest"
};

static std::string pcommtestclasses[] = { "P_MPICommunicatorTest" };

class CommunicatorTest
{
public:

    CommunicatorTest( const scai::dmemo::communicator::CommunicatorKind ct );
    ~CommunicatorTest();

    void CommunicatorCtrTest();

    template<typename ValueType> void swapTest() __attribute__( ( noinline ) );
    template<typename ValueType> void gatherVTest() __attribute__( ( noinline ) );
    template<typename ValueType> void gatherTest() __attribute__( ( noinline ) );
    template<typename ValueType> void scatterVTest() __attribute__( ( noinline ) );
    template<typename ValueType> void scatterTest() __attribute__( ( noinline ) );
    template<typename ValueType> void bcastTest() __attribute__( ( noinline ) );
    template<typename ValueType> void shiftASyncTest() __attribute__( ( noinline ) );
    template<typename ValueType> void shiftTest() __attribute__( ( noinline ) );
    template<typename ValueType> void updateHaloTest() __attribute__( ( noinline ) );

    void bcastStringTest();
    void buildHaloTest();
    void allocatePlanTest();
    void computeOwnersTest();
    void writeAtTest();

    void runTests();

private:
    scai::dmemo::communicator::CommunicatorKind mCommunicatorType;

    scai::dmemo::CommunicatorPtr comm;
    PartitionId rank;
    PartitionId size;
};

#define COMMUNICATORTEST_COMMONTESTCASES( testinstance )                                     \
    {   COMMONTESTCASEINVOKER_TEMPLATE( testinstance, swapTest, float );                     \
        COMMONTESTCASEINVOKER_TEMPLATE( testinstance, gatherTest, float );                   \
        COMMONTESTCASEINVOKER_TEMPLATE( testinstance, gatherVTest, float );                  \
        COMMONTESTCASEINVOKER_TEMPLATE( testinstance, scatterTest, float );                  \
        COMMONTESTCASEINVOKER_TEMPLATE( testinstance, scatterVTest, float );                 \
        COMMONTESTCASEINVOKER_TEMPLATE( testinstance, bcastTest, float );                    \
        COMMONTESTCASEINVOKER_TEMPLATE( testinstance, shiftASyncTest, float );               \
        COMMONTESTCASEINVOKER_TEMPLATE( testinstance, shiftTest, float );                    \
        COMMONTESTCASEINVOKER_TEMPLATE( testinstance, updateHaloTest, float );               \
        COMMONTESTCASEINVOKER_TEMPLATE( testinstance, swapTest, double );                    \
        COMMONTESTCASEINVOKER_TEMPLATE( testinstance, gatherTest, double );                  \
        COMMONTESTCASEINVOKER_TEMPLATE( testinstance, gatherVTest, double );                 \
        COMMONTESTCASEINVOKER_TEMPLATE( testinstance, scatterTest, double );                 \
        COMMONTESTCASEINVOKER_TEMPLATE( testinstance, scatterVTest, double );                \
        COMMONTESTCASEINVOKER_TEMPLATE( testinstance, bcastTest, double );                   \
        COMMONTESTCASEINVOKER_TEMPLATE( testinstance, shiftASyncTest, double );              \
        COMMONTESTCASEINVOKER_TEMPLATE( testinstance, shiftTest, double );                   \
        COMMONTESTCASEINVOKER_TEMPLATE( testinstance, updateHaloTest, double );              \
        COMMONTESTCASEINVOKER( testinstance, bcastStringTest );                              \
        COMMONTESTCASEINVOKER( testinstance, buildHaloTest );                                \
        COMMONTESTCASEINVOKER( testinstance, allocatePlanTest );                             \
        COMMONTESTCASEINVOKER( testinstance, computeOwnersTest );                            \
        COMMONTESTCASEINVOKER( testinstance, CommunicatorCtrTest ); }
