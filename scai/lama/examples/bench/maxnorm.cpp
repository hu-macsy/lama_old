/**
 * @file maxnorm.cpp
 *
 * @license
 * Copyright (c) 2009-2018
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the SCAI framework LAMA.
 *
 * LAMA is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * LAMA is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
 * more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with LAMA. If not, see <http://www.gnu.org/licenses/>.
 * @endlicense
 *
 * @brief Test of maxnorm for all valuetypes
 * @author Eric Schricker
 * @date 21.03.2016
 */

#include <iostream>
#include <iomanip>

#include <scai/lama.hpp>

// _Matrix & vector related includes
#include <scai/lama/matrix/all.hpp>
#include <scai/lama/DenseVector.hpp>

#include <scai/dmemo/Distribution.hpp>
#include <scai/dmemo/BlockDistribution.hpp>

#include <scai/common/Walltime.hpp>

using namespace scai::lama;
using namespace scai::hmemo;
using namespace scai::dmemo;
using namespace std;
using scai::IndexType;
using scai::common::Walltime;

template<typename ValueType>
static void bench( IndexType size )
{
    auto dist = std::make_shared<BlockDistribution>( size );
    DenseVector<ValueType> x( dist, 7 );
    double tmpTime = Walltime::get();
    x.maxNorm();
    tmpTime = Walltime::get() - tmpTime;
    std::cout << std::left << std::setw( 15 ) << std::setfill( ' ' ) << std::setprecision( 3 ) << tmpTime;
}

int main()
{
    IndexType sizes[] = { 128, 256, 512, 1024, 2048, 4096, 8192,
                          16384, 32768, 65536, 131072, 262144,
                          524288, 1048576, 2097152, 4194304,
                          8388608, 16777216, 33554432
                        };
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

    for ( IndexType i = 0; i < n; ++i )
    {
        std::cout << std::left << std::setw( 15 ) << std::setfill( ' ' ) << std::setprecision( 3 ) << sizes[i];
#define DO_BENCH( ValueType ) bench<ValueType>( sizes[i] );
        // do the benchmark for each supported A type
        SCAI_COMMON_LOOP( DO_BENCH, SCAI_NUMERIC_TYPES_HOST )
#undef DO_BENCH
        std::cout << std::endl;
    }
}
