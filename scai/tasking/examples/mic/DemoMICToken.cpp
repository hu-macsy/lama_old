/**
 * @file tasking/examples/mic/DemoMICToken.cpp
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
 * @brief Example of using a MIC sync token
 *        /
 * @author Thomas Brandes
 * @date 30.03.2016
 */

#include <scai/logging.hpp>
#include <scai/tasking/mic/MICSyncToken.hpp>
#include <scai/common/Walltime.hpp>

#include <iostream>
#include <unistd.h>
#include <cstdlib>

using namespace scai;
using namespace scai::common;
using namespace scai::tasking;

const int NSIZE = 700;

__declspec( target( mic, host ) )
extern void compute( size_t NSIZE );

void start_mic_compute( MICSyncToken& token )
{
    int signal_value = 1;
    int device = token.getDevice();
    int n = NSIZE;
#pragma offload target( mic : device ), in( n ), signal(signal_value)
    {
        compute( n );
    }
    token.setSignal( signal_value );
}

int main()
{
    double t;
    int device = 0;
    int n = NSIZE;
    // measure coprocessor
    t = Walltime::get();
#pragma offload target(mic : device)
    {
        compute( n );
    }
    t = Walltime::get() - t;
    std::cout << "MIC run (1st): " << t << " seconds" << std::endl;
    // measure coprocessor again, no init device
    t = Walltime::get();
#pragma offload target(mic : device)
    {
        compute( n );
    }
    t = Walltime::get() - t;
    std::cout << "MIC run (2nd): " << t << " seconds" << std::endl;
    t = Walltime::get();
    compute( n );
    t = Walltime::get() - t;
    std::cout << "Host activity = " << t << " seconds" << std::endl;
    // LAMA approach using SyncToken
    t = Walltime::get();
    {
        MICSyncToken token( device );
        start_mic_compute( token );
        compute( n );
        // implicit wait at end of this scope
    }
    t = Walltime::get() - t;
    std::cout << "Host + Device activity = " << t << " seconds" << std::endl;
}
