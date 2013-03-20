/**
 * @file DistributionTest.cpp
 *
 * @license
 * Copyright (c) 2011
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
 * @author Alexander Büchel, Thomas Brandes
 * @date 30.07.2012
 * $Id$
 */

#include <test/distributed/DistributionTest.hpp>

#include <test/TestMacros.hpp>

LAMA_COMMON_TEST_CASE( DistributionTest, localSizeTest );
{
    // general test that must work for all distributions
    // the local sizes summed up over all partitions must be the global size

    IndexType sumLocalSizes = mDistributionPtr->getCommunicator().sum( mDistributionPtr->getLocalSize() );
    BOOST_CHECK_EQUAL( mDistributionPtr->getGlobalSize(), sumLocalSizes );
}
LAMA_COMMON_TEST_CASE_END()
;

LAMA_COMMON_TEST_CASE( DistributionTest, local2GlobalTest );
{
    for ( IndexType i = 0; i < mDistributionPtr->getGlobalSize(); i++ )
    {
        if ( mDistributionPtr->isLocal( i ) )
        {
            BOOST_CHECK_EQUAL( i, mDistributionPtr->local2global( mDistributionPtr->global2local(i) ) );
        }
        else
        {
            BOOST_CHECK_EQUAL( nIndex, mDistributionPtr->global2local(i) );
        }
    }
}
LAMA_COMMON_TEST_CASE_END()
;

LAMA_COMMON_TEST_CASE( DistributionTest, global2LocalTest );
{
    for ( IndexType i = 0; i < mDistributionPtr->getLocalSize(); i++ )
    {
        BOOST_CHECK_EQUAL( i, mDistributionPtr->global2local( mDistributionPtr->local2global(i) ) );
    }
}
LAMA_COMMON_TEST_CASE_END()
;

LAMA_COMMON_TEST_CASE( DistributionTest, writeAtTest );
{
    LAMA_WRITEAT_PTR_TEST( mDistributionPtr );
}
LAMA_COMMON_TEST_CASE_END()
;

LAMA_COMMON_TEST_CASE( DistributionTest, printDistributionVector );
{
    //todo: does not test the content of these files
    std::string s;
    s = "distribution";
    mDistributionPtr->printDistributionVector( s );
    std::remove( ( s + ".part" ).c_str() );
}
LAMA_COMMON_TEST_CASE_END()
;

LAMA_COMMON_TEST_CASE_RUNNER( DistributionTest )
{
    writeAtTest();
    localSizeTest();
    local2GlobalTest();
    global2LocalTest();
    printDistributionVector();
}
