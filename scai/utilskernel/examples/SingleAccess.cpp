/**
 * @file SingleAccess.cpp
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
 *
 * Other Usage
 * Alternatively, this file may be used in accordance with the terms and
 * conditions contained in a signed written agreement between you and
 * Fraunhofer SCAI. Please contact our distributor via info[at]scapos.com.
 * @endlicense
 *
 * @brief Benchmarking of single remote read and write operations
 * @author Thomas Brandes
 * @date 03.07.2015
 */

#include <scai/utilskernel.hpp>

#include <scai/common/Walltime.hpp>

using namespace std;
using namespace scai;
using namespace hmemo;
using namespace utilskernel;

int main()
{
    typedef RealType ValueType;

    if ( ! Context::canCreate( Context::CUDA ) )
    {
        cout << "Example program skipped, no CUDA available." << endl;
    }

    ContextPtr host = Context::getHostPtr();
    ContextPtr gpu  = Context::getContextPtr( Context::CUDA );

    const IndexType Nh = 50;
    const IndexType N = 2 * Nh;

    LArray<ValueType> hostA( N, 5, host );
    LArray<ValueType> gpuA( N, 2, gpu );

    for ( IndexType i = 0; i < N; i+=2 )
    {
        gpuA[i] = hostA[i];
    }

    IndexType sum = gpuA.sum();

    cout << "Sum = " << sum << " should be " << 2 * Nh + 5 * Nh  << endl;

    {
        // invalidate copy of gpuA on host

        WriteAccess<ValueType> write( gpuA, gpu );
    }

    for ( IndexType i = 1; i < N; i+=2 )
    {
        hostA[i] = gpuA[i];
    }

    sum = hostA.sum();

    cout << "Sum = " << sum << " should be " << 2 * Nh + 5 * Nh  << endl;
}
