/**
 * @file RemoteAccess.cpp
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

#include <scai/hmemo.hpp>

#include <scai/common/Walltime.hpp>

using namespace std;
using namespace scai;
using namespace hmemo;
using common::ContextType;

int main()
{
    // Note: in hmemo we can use any value type for HArray

    typedef double ValueType;

    if ( ! Context::canCreate( ContextType::CUDA ) )
    {
        cout << "Example program skipped, no CUDA available." << endl;
        return 0;
    }

    ContextPtr host = Context::getHostPtr();
    ContextPtr gpu  = Context::getContextPtr( ContextType::CUDA );

    const IndexType Nh = 50;
    const IndexType N = 2 * Nh;

    HArray<ValueType> hostA( N, 5, host );
    HArray<ValueType> gpuA( N, 2, gpu );

    for ( IndexType i = 0; i < N; i += 2 )
    {
        ReadAccess<ValueType> readA( hostA, host );
        WriteAccess<ValueType> writeA( gpuA, gpu );
        ValueType elem = readA[i];
        writeA.setValue( elem, i );
    }

    ValueType sum = 0;

    {
        ReadAccess<double> readA( gpuA, host );

        for ( IndexType i = 0; i < N; ++i )
        {
            sum += readA[i];
        }
    }

    cout << "Sum = " << sum << " should be " << 2 * Nh + 5 * Nh  << endl;

    {
        // invalidate copy of gpuA on host

        WriteAccess<ValueType> write( gpuA, gpu );
    }

    for ( IndexType i = 1; i < N; i += 2 )
    {
        ReadAccess<ValueType> readA( gpuA, gpu );
        WriteAccess<ValueType> writeA( hostA, host );
        ValueType elem;
        readA.getValue( elem, i );
        cout << "val @ " << i << " = " << elem << endl;
        writeA[i] = elem;
    }

    sum = 0;

    {
        ReadAccess<double> readA( hostA, host );

        for ( IndexType i = 0; i < N; ++i )
        {
            sum += readA[i];
        }
    }

    cout << "Sum = " << sum << " should be " << 2 * Nh + 5 * Nh  << endl;
}
