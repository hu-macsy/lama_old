/**
 * @file examples/cuda/DemoCUDACtx.cpp
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
 * @brief ToDo: Missing description in ./examples/cuda/DemoCUDACtx.cpp
 * @author Thomas Brandes
 * @date 09.03.2016
 */

#include <scai/common/cuda/CUDACtx.hpp>
#include <scai/common/cuda/CUDAError.hpp>

#include <scai/common/Settings.hpp>

#include <iostream>

using namespace scai;
using namespace common;

int main( int argc, const char** argv )
{
    // at least --SCAI_DEVICE=id may be specified
    Settings::parseArgs( argc, argv );
    int nr = 0;   // take this as default
    Settings::getEnvironment( nr, "SCAI_DEVICE" );
    CUDACtx context( nr );
    // Note: no context access required for queries of the device
    char deviceName[256];
    SCAI_CUDA_DRV_CALL( cuDeviceGetName( deviceName, 256, context.getCUdevice() ), "cuDeviceGetName" );
    std::cout << "CUDACtx( device = " << context.getDeviceNr() << " ), name of device = " << deviceName << std::endl;
}
