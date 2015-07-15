/**
 * @file Devices.cu   
 *
 * @license
 * Copyright (c) 2009-2015
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 * @endlicense
 *
 * @brief Example program to show how to query and access CUDA devices.
 * @author: Thomas Brandes
 * @date 15.06.2015
 **/

#include <iostream>

#include <memory/memory.hpp>

#include <common/Exception.hpp>

#include <logging/logging.hpp>

using namespace std;
using namespace memory;

void sub( ContextPtr cudaContext )
{
    try
    {
        MemoryPtr mem = cudaContext->getMemory();
        cout << "CUDA context " << *cudaContext << " has mem = " << *mem << endl;
    }
    catch ( common::Exception& ex )
    {
        cout << "CUDA memory for " << *cudaContext << " failed" << endl;
    }

    try
    {
        MemoryPtr mem = cudaContext->getHostMemory();
        cout << "CUDA context " << *cudaContext << " has host mem = " << *mem << endl;
    }
    catch ( common::Exception& ex )
    {
        cout << "CUDA host memory for " << *cudaContext << " failed" << endl;
    }

}

int main()
{
    for ( int deviceNr = 0; deviceNr < 8; ++ deviceNr )
    {
        cout << "try to get " << context::CUDA << " context from factory" << endl;

        try 
        {
            ContextPtr cudaContext = Context::create( context::CUDA, deviceNr );
            cout << "cudaContext for device " << deviceNr << " = " << *cudaContext << endl;
            sub( cudaContext );
        }
        catch ( common::Exception& ex )
        {
            cout << "CUDA device " << deviceNr << " is not available" << endl;
        }
    }
}

