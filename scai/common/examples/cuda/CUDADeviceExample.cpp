
#include <scai/common/cuda/CUDADevice.hpp>
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
   
    CUDADevice device( nr );

    // Note: no context access required for queries of the device 

    char deviceName[256];

    SCAI_CUDA_DRV_CALL( cuDeviceGetName( deviceName, 256, device.getCUdevice() ), "cuDeviceGetName" );

    std::cout << "Name of device = " << deviceName << std::endl;
}
