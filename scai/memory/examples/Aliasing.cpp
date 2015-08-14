/**
 * @file Aliasing.cpp
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
 * @brief Demo on aliasing problem with LAMA arrays
 * @author: Thomas Brandes
 * @date 03.07.2015
 **/

#include <scai/memory/Context.hpp>
#include <scai/memory/ReadAccess.hpp>
#include <scai/memory/WriteOnlyAccess.hpp>
#include <scai/common/Exception.hpp>

using namespace memory;

/* --------------------------------------------------------------------- */

SCAI_LOG_DEF_LOGGER( logger, "ContextTest" )

using namespace memory;
using common::IndexType;

typedef LAMAArray<double> Array;

void add ( Array& res, const Array& a, const Array& b )
{
    COMMON_ASSERT_LE( res.size(), a.size(), "size mismatch" )
    COMMON_ASSERT_LE( res.size(), b.size(), "size mismatch" )

    IndexType n = res.size();

    ContextPtr host = Context::getContextPtr( context::Host );

    WriteOnlyAccess<double> write( res, host );
    ReadAccess<double>read1( a, host );
    ReadAccess<double>read2( b, host );
 
    double* resPtr = write.get();
    const double* aPtr = read1.get();
    const double* bPtr = read2.get();

    for ( IndexType i = 0; i < n; ++i )
    {
        resPtr[i] = aPtr[i] + bPtr[i];
    }
}

int main()
{
    Array a( 10 );
    Array b( 10 , 1.0 );
    Array c( 10 , 1.0 );

    add( a, b, c ); // this is okay
    add( a, a, b ); // this crashed in earlier versions
}

