/**
 * @file spy.cpp
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
 * @brief Example program to spy sparse structure of a CSR matrix.
 * @author Thomas Brandes
 * @date 17.10.2013
 * @since 1.0.0
 */

#include "Bitmap.hpp"

#include <lama/storage/CSRStorage.hpp>
#include <lama/HostReadAccess.hpp>

using namespace lama;

int main( int argc, char** argv )
{
    CSRStorage<double> matrix;

    if ( argc < 2 )
    {
        std::cerr << "Missing filename for input matrix" << std::endl;
        std::cerr << "spy matrix_filename [ width [ height [ scale ] ] ]" << std::endl;
        exit(1);
    }

    const char* filename = argv[1];

    int nRows = 800;

    if ( argc > 2 )
    {
        sscanf( argv[2], "%d",  &nRows );
    }

    int nColumns = nRows;

    if ( argc > 3 )
    {
        sscanf( argv[3], "%d",  &nColumns );
    }

    int nZoom = 1;

    if ( argc > 4 )
    {
        sscanf( argv[4], "%d",  &nZoom );
    }

    matrix.readFromFile( filename );

    const LAMAArray<IndexType>& ia = matrix.getIA();
    const LAMAArray<IndexType>& ja = matrix.getJA();
    const LAMAArray<double>& values = matrix.getValues();

    HostReadAccess<IndexType> csrIA( ia );
    HostReadAccess<IndexType> csrJA( ja );
    HostReadAccess<double> csrValues( values );

    std::cout << "Write png of size " << nRows << " x " << nColumns << ", zoom = " << nZoom << std::endl;

    Bitmap pic( nRows, nColumns, nZoom );

    pic.setColor( 240, 120, 0 );  // color for smallest value
    // pic.setColor( 0, 0, 255 );    // color for largetst value

    pic.drawCSR( matrix.getNumRows(), matrix.getNumColumns(), csrIA.get(), csrJA.get(), csrValues.get() );

    const std::string out_filename = "lama.png";

    pic.write_png_file( out_filename.c_str() );

    std::cout << "png files has been written as " << out_filename << std::endl;
}
