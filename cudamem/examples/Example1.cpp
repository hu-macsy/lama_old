
#include <memory/LAMAArray.hpp>
#include <memory/Context.hpp>

#include <iostream>

using namespace memory;

int main()
{
    ContextPtr cudaContext = Context::getContext( context::CUDA );

    std::cout << "cudaContext = " << *cudaContext << std::endl;

    LAMAArray<double> data( 100, 1.0 );
    
    std::cout << "data = " << data << std::endl;
}

