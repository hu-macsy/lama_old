/**
 * @file NormTest.hpp
 *
 * @license
 * Copyright (c) 2009-2016
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the Library of Accelerated Math Applications (LAMA).
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
 * @brief Contains the implementation of the class NormTest
 * @author Alexander Büchel, Micha
 * @date 03.02.2012
 */

#include <boost/test/unit_test.hpp>

#include <scai/lama/norm/Norm.hpp>
#include <scai/lama/test/TestMacros.hpp>

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

    NormTest( const scai::lama::Norm& norm )
        : mNorm( norm )
    {
    }

    void positiveHomogeneityTest();
    void triangleInequalityTest();
    void ZeroVectorTest();
    void runTests();

    const scai::lama::Norm& mNorm;
};

#define NORMTEST_COMMONTESTCASES( testinstance )                        \
    {   COMMONTESTCASEINVOKER( testinstance, positiveHomogeneityTest );     \
        COMMONTESTCASEINVOKER( testinstance, triangleInequalityTest );      \
        COMMONTESTCASEINVOKER( testinstance, ZeroVectorTest ); }
