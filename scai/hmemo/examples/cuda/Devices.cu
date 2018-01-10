/**
 * @file hmemo/examples/cuda/Devices.cu
 *
 * @license
 * Copyright (c) 2009-2017
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
 * @brief Example program to show how to query and access CUDA devices.
 * @author Thomas Brandes
 * @date 15.06.2015
 */

#include <iostream>

#include <scai/hmemo.hpp>

#include <scai/common/macros/throw.hpp>

#include <scai/logging.hpp>

using namespace std;
using namespace scai;
using namespace scai::hmemo;

void sub( ContextPtr cudaContext )
{
    try
    {
        MemoryPtr mem = cudaContext->getMemoryPtr();
        cout << "CUDA context " << *cudaContext << " has mem = " << *mem << endl;
    }
    catch ( common::Exception& ex )
    {
        cout << "CUDA memory for " << *cudaContext << " failed" << endl;
    }

    try
    {
        MemoryPtr mem = cudaContext->getHostMemoryPtr();
        cout << "CUDA context " << *cudaContext << " has host mem = " << *mem << endl;
    }
    catch ( common::Exception& ex )
    {
        cout << "CUDA host memory for " << *cudaContext << " failed" << endl;
    }
}

int main()
{
    if ( Context::canCreate( common::ContextType::CUDA ) )
    {
        cout << "Factory can create CUDA context, registered" << endl;
    }
    else
    {
        cout << "Factory cannot create CUDA context, not registed" << endl;
        exit( -1 );  // continuation makes no sense
    }

    for ( int deviceNr = 0; deviceNr < 8; ++ deviceNr )
    {
        cout << "try to get " << common::ContextType::CUDA << " context from factory" << endl;

        try
        {
            ContextPtr cudaContext = Context::create( common::ContextType::CUDA, deviceNr );
            cout << "cudaContext for device " << deviceNr << " = " << *cudaContext << endl;
            sub( cudaContext );
        }
        catch ( common::Exception& ex )
        {
            cout << "CUDA device " << deviceNr << " is not available" << endl;
        }
    }
}

