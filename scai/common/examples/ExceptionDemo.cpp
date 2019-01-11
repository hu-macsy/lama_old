/**
 * @file common/examples/ExceptionDemo.cpp
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
 * @brief Demo example program using Exception and ASSERT macros.
 * @author Thomas Brandes
 * @date 19.06.2015
 */

#include <scai/common/macros/throw.hpp>
#include <scai/common/macros/assert.hpp>

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

