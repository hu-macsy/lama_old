/**
 * @file common/examples/ExceptionDemo.cpp
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
 * @brief Demo example program using Exception and ASSERT macros.
 *
 * @author Thomas Brandes
 * @date 19.06.2015
 */

#include <scai/common/Exception.hpp>

#include <iostream>

void sub( int val )
{
    if ( val < 0 ) 
    {
        COMMON_THROWEXCEPTION( "sub: val must not be negative, val = " << val )
    }

    SCAI_ASSERT ( val % 2 == 0, "val = " << val << " must be even" )

    SCAI_ASSERT_LT( val, 10, "val = " << val << " must be less than 10" )

    SCAI_ASSERT_EQUAL( val, 4, "None" )
}

int main()
{
    int vals[] = { -1, 5, 14, 6, 4 };
 
    int nargs = sizeof( vals ) / sizeof( int );

    for ( int i = 0; i < nargs; ++ i )
    {
        try
        {
            sub( vals[i] );
            std::cout << "Call of sub( " << vals[i] << ") terminated correctly" << std::endl;
        }
        catch ( const std::exception& exception )
        {
            // Importation: exception is a reference, so we get the routine of common::Exception
    
            std::cout << "Got exception: " << exception.what() << std::endl;
        }
    }
}

