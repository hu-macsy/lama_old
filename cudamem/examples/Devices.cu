
#include <iostream>

#include <common/Exception.hpp>

#include <logging/logging.hpp>

// using namespace std;
// using namespace memory;

int main()
{
    for ( int deviceNr = 0; deviceNr < 8; ++ deviceNr )
    {
        cout << "try to get " << context::CUDA << " context from factory" << endl;

        try 
        {
            ContextPtr cudaContext = Context::create( context::CUDA, deviceNr );
            cout << "cudaContext for device " << deviceNr << " = " << *cudaContext << endl;
        }
        catch ( common::Exception& ex )
        {
            cout << "CUDA device " << deviceNr << " is not available" << endl;
        }
    }
}

