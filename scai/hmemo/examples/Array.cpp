/**
 * @file hmemo/examples/Array.cpp
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
 * @brief Contains the implementation of the class HArrayTest.
 * @author Thomas Brandes, Lauretta Schubert
 * @date 18.04.2012
 */

#include <scai/hmemo/HArray.hpp>
#include <scai/hmemo/WriteAccess.hpp>
#include <scai/hmemo/ReadAccess.hpp>

#include <scai/logging.hpp>

#include <iostream>

using namespace std;
using namespace scai;
using namespace scai::hmemo;

SCAI_LOG_DEF_LOGGER( logger, "MemoryTest" )

template<typename T>
void sumArray( const HArray<T>& array )
{
    SCAI_LOG_INFO( logger, "read access on " << array );
    ContextPtr contextPtr = Context::getHostPtr();
    ReadAccess<T> readAccess( array, contextPtr );
    const T* data = readAccess.get();
    T sum = 0;

    for ( IndexType i = 0; i < array.size(); ++i )
    {
        sum += data[i];
    }

    cout << "Sum of " << array << " = " << sum << endl;
}

template<typename T>
void writeArray( HArray<T>& array )
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

struct SSS
{
    int X;
    double Y;

    SSS( int val )
    {
        X = val;
        Y = 0;
    }

    SSS( int valX, double valY )
    {
        X = valX;
        Y = valY;
    }

    SSS& operator +=( const SSS& other )
    {
        X += other.X;
        Y += other.Y;
        return *this;
    }

};

std::ostream& operator<<( std::ostream& stream, const SSS& object )
{
    stream << "( " << object.X << ", " << object.Y << " )";
    return stream;
}

int main()
{
    SCAI_LOG_THREAD( "Main" )
    HArray<IndexType> lamaArray; // default, not allocated at all
    sumArray( lamaArray );
    writeArray( lamaArray );
    sumArray( lamaArray );
    SSS val( 1, 0.5 );
    val.X = 1;
    val.Y = 0.5;
    HArray<SSS> sssArray( 10, val );
    sumArray( sssArray );
}
