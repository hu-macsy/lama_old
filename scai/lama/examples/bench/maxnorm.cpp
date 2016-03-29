/**
 * @file maxnorm.cpp
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
 * @brief Test of maxnorm for all valuetypes
 * @author Eric Schricker
 * @date 21.03.2016
 * @since 2.0.0
 */

#include <iostream>
#include <iomanip>

#include <scai/lama.hpp>

// Matrix & vector related includes
#include <scai/lama/expression/all.hpp>
#include <scai/lama/matrix/all.hpp>
#include <scai/lama/DenseVector.hpp>

#include <scai/dmemo/Distribution.hpp>
#include <scai/dmemo/BlockDistribution.hpp>

#include <scai/common/Walltime.hpp>

using namespace scai::lama;
using namespace scai::hmemo;
using namespace scai::dmemo;
using namespace std;
using scai::common::Walltime;

template<typename ValueType>
static void bench( IndexType size )
{
    CommunicatorPtr comm = Communicator::getCommunicatorPtr();
    DistributionPtr dist( new BlockDistribution( size, comm ));

    DenseVector<ValueType> x( dist );

    x = 7.0;

    double tmpTime = Walltime::get();
    x.maxNorm();
    tmpTime = Walltime::get() - tmpTime;
    std::cout << std::left << std::setw( 15 ) << std::setfill( ' ' ) << std::setprecision(3) << tmpTime;
}

int main()
{
    IndexType sizes[] = { 128, 256, 512, 1024, 2048, 4096, 8192,
                          16384, 32768, 65536, 131072, 262144,
                          524288, 1048576, 2097152, 4194304,
                          8388608, 16777216, 33554432 };
    IndexType n = sizeof( sizes ) / sizeof( IndexType );

    std::cout << std::left << std::setw( 15 ) << std::setfill( ' ' ) << "Size";
    std::cout << std::left << std::setw( 15 ) << std::setfill( ' ' ) << "float";
    std::cout << std::left << std::setw( 15 ) << std::setfill( ' ' ) << "double";
    std::cout << std::left << std::setw( 15 ) << std::setfill( ' ' ) << "long double";
    std::cout << std::left << std::setw( 15 ) << std::setfill( ' ' ) << "ComplexFloat";
    std::cout << std::left << std::setw( 15 ) << std::setfill( ' ' ) << "ComplexDouble";
    std::cout << std::left << std::setw( 15 ) << std::setfill( ' ' ) << "ComplexLongDouble";
    std::cout << std::endl;

    std::cout << std::setw( 110 ) << std::setfill( '-' ) << "-";
    std::cout << std::endl;

    for( IndexType i = 0; i < n; ++i )
    {

        std::cout << std::left << std::setw( 15 ) << std::setfill( ' ' ) << std::setprecision(3) << sizes[i];
        bench<float>( sizes[i] );
        bench<double>( sizes[i] );
        bench<long double>( sizes[i] );
        bench<ComplexFloat>( sizes[i] );
        bench<ComplexDouble>( sizes[i] );
        bench<ComplexLongDouble>( sizes[i] );
        std::cout << std::endl;
    }

}

