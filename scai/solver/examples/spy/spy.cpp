/**
 * @file spy.cpp
 *
 * @license
 * Copyright (c) 2009-2017
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the SCAI framework LAMA.
 *
 * LAMA is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Affero General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * LAMA is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for
 * more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with LAMA. If not, see <http://www.gnu.org/licenses/>.
 *
 * Other Usage
 * Alternatively, this file may be used in accordance with the terms and
 * conditions contained in a signed written agreement between you and
 * Fraunhofer SCAI. Please contact our distributor via info[at]scapos.com.
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

typedef RealType ValueType;

int main( int argc, char** argv )
{
    CSRStorage<ValueType> matrix;

    if ( argc < 2 )
    {
        std::cerr << "Missing filename for input matrix" << std::endl;
        std::cerr << "spy matrix_filename [ width [ height [ scale ] ] ]" << std::endl;
        exit( 1 );
    }

    const char* filename = argv[1];

    matrix.readFromFile( filename );

    IndexType height = matrix.getNumRows();
    IndexType width = matrix.getNumColumns();

    while ( width > 2048 || height > 2048 )
    {
        width = width / 2;
        height = height / 2;
    }

    if ( argc > 2 )
    {
        sscanf( argv[2], "%d",  &width );
        height = width;
    }

    if ( argc > 3 )
    {
        sscanf( argv[3], "%d",  &height );
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
    const std::string out_filename = "lama.png";
    pic.write( out_filename );
    std::cout << "png files has been written as " << out_filename << std::endl;
}
