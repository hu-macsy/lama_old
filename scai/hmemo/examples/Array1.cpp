/**
 * @file Array1.cpp
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
 * @brief Contains the implementation of the class LAMAArrayTest.
 * @author: Thomas Brandes, Lauretta Schubert
 * @date 18.04.2012
 **/

#include <scai/hmemo/LAMAArray.hpp>
#include <scai/hmemo/WriteAccess.hpp>
#include <scai/hmemo/ReadAccess.hpp>

#include <scai/logging.hpp>

#include <iostream> 

using namespace std;
using namespace scai;
using namespace scai::hmemo;

SCAI_LOG_DEF_LOGGER( logger, "MemoryTest" )

template<typename T>
void sumArray( const LAMAArray<T>& array )
{
    SCAI_LOG_INFO( logger, "read access on " << array );

    ContextPtr contextPtr = Context::getHostPtr();

    ReadAccess<T> readAccess( array, contextPtr );
   
    const T* data = readAccess.get();
 
    T sum = 0;

    for ( int i = 0; i < array.size(); ++i )
    {
        sum += data[i];
    }

    cout << "Sum of " << array << " = " << sum << endl;
}

template<typename T>
void writeArray( LAMAArray<T>& array )
{
    SCAI_LOG_INFO( logger, "make write test access on empty array\n" );

    ContextPtr contextPtr = Context::getHostPtr();

    WriteAccess<T> writeAccess( array, contextPtr );

    writeAccess.resize( 10 );

    T* data = writeAccess.get();

    // data is on host, so we can work directly on it

    for ( IndexType i = 0; i < 10; i++ )
    {
        data[i] = 3;
    }
}

struct SSS {
int X; double Y;
};

int main()
{
    SCAI_LOG_THREAD( "Main" )

    LAMAArray<IndexType> lamaArray; // default, not allocated at all

    sumArray( lamaArray );
    writeArray( lamaArray );
    sumArray( lamaArray );

    LAMAArray<SSS> sssArray( 10 );
}
