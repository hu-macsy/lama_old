#include <CL/opencl.h>
#include <iostream>

extern "C" cl_int oclGetPlatformID(cl_platform_id* clSelectedPlatformID);

////////////////////////////////////////////////////////////////////////////////
// Program main
////////////////////////////////////////////////////////////////////////////////

using namespace std;

int main( int argc, char** argv)
{
    cl_uint num_platforms;

    cl_int ciErrNum = clGetPlatformIDs ( 0, NULL, &num_platforms );

    if ( ciErrNum != CL_SUCCESS )
    {
        cerr << "clGetPlatformIDs failed" << endl;
        return -1;
    }

    cout << "#OpenCL platforms = " << num_platforms << endl;

}
