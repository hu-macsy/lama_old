/**
 * @file DemoDistribution.cpp
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
 * @brief Demo program for distribution + communication
 * @author: Thomas Brandes
 * @date 10.02.2016
 **/

#include <scai/dmemo.hpp>

using namespace scai::dmemo;

int main()
{
    SCAI_LOG_THREAD( "Main" )

    // get the default communicator (usually MPI if it has been enabled, or set by SCAI_COMMUNICATOR

    CommunicatorPtr comm = Communicator::getCommunicator();

    IndexType size = 71;

    float weight = 1.0;

    DistributionPtr dist ( Distribution::getDistribution( "CYCLIC", comm, size, weight ) );

    // Note: distribution pointers are always const pointers, so distributions can never be changed

    std::cout << *comm << ", dist = " << *dist << std::endl;
}
