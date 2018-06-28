/**
 * @file SingleAccess.cpp
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
    typedef DefaultReal ValueType;

    if ( ! Context::canCreate( common::ContextType::CUDA ) )
    {
        cout << "Example program skipped, no CUDA available." << endl;
        return 0;
    }

    ContextPtr host = Context::getHostPtr();
    ContextPtr gpu  = Context::getContextPtr( common::ContextType::CUDA );

    const IndexType Nh = 50;
    const IndexType N = 2 * Nh;

    HArray<ValueType> hostA( N, 5, host );
    HArray<ValueType> gpuA( N, 2, gpu );

    for ( IndexType i = 0; i < N; i += 2 )
    {
        gpuA[i] = hostA[i];
    }

    ValueType sum = HArrayUtils::sum( gpuA );

    cout << "Sum = " << sum << " should be " << 2 * Nh + 5 * Nh  << endl;

    {
        // invalidate copy of gpuA on host

        WriteAccess<ValueType> write( gpuA, gpu );
    }

    for ( IndexType i = 1; i < N; i += 2 )
    {
        hostA[i] = gpuA[i];
    }

    sum = HArrayUtils::sum( hostA );

    cout << "Sum = " << sum << " should be " << 2 * Nh + 5 * Nh  << endl;
}
