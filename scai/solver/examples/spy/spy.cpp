/**
 * @file spy.cpp
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
 * @brief Example program to spy sparse structure of a CSR matrix.
 * @author Thomas Brandes
 * @date 17.10.2013
 */

#include "Bitmap.hpp"

#include <scai/lama/storage/CSRStorage.hpp>
#include <scai/hmemo/ReadAccess.hpp>

using namespace scai::lama;
using namespace scai::hmemo;

typedef DefaultReal ValueType;

int main( int argc, char** argv )
{
    CSRStorage<ValueType> matrix;

    if ( argc < 3 )
    {
        std::cerr << "Missing filename for input matrix" << std::endl;
        std::cerr << "spy matrix_filename bitmap_filename [ width [ height ] ]" << std::endl;
        exit( 1 );
    }

    std::string matrixFileName( argv[1] );
    std::string imageFileName( argv[2] );

    matrix.readFromFile( matrixFileName );

    IndexType height = matrix.getNumRows();
    IndexType width = matrix.getNumColumns();

    if ( argc > 3 )
    {
        std::istringstream input( argv[3] );
        input >> width;
        height = width;
    
        if ( argc > 4 )
        {
            std::istringstream input( argv[4] );
            input >> height;
        }

        SCAI_ASSERT_GE_ERROR( matrix.getNumRows(), width, "width cannot be greater than #cols in matrix" );
        SCAI_ASSERT_GE_ERROR( matrix.getNumColumns(), height, "height cannot be greater than #rows in matrix" );
    }

    else
    {
        height = matrix.getNumRows();
        width = matrix.getNumColumns();

        while ( width > 2048 || height > 2048 )
        {
            width = width / 2;
            height = height / 2;
        }
    }

    const HArray<IndexType>& ia = matrix.getIA();
    const HArray<IndexType>& ja = matrix.getJA();
    const HArray<ValueType>& values = matrix.getValues();

    ReadAccess<IndexType> csrIA( ia );
    ReadAccess<IndexType> csrJA( ja );
    ReadAccess<ValueType> csrValues( values );

    std::cout << "Write png of size " << height << " x " << width << std::endl;

    Bitmap pic( height, width );
    pic.drawCSR( matrix.getNumRows(), matrix.getNumColumns(), csrIA.get(), csrJA.get(), csrValues.get() );
    pic.write( imageFileName );

    std::cout << "Done: written matrix " << matrixFileName 
              << " into image file " << imageFileName 
              << ", " << width << " x " << height << std::endl;
}
