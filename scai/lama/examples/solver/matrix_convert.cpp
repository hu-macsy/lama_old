/**
 * @file matrix_convertor.cpp
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
 * @brief Conversion program between binary and matrix market
 * @author Thomas Brandes
 * @date 20.12.2015
 */

// Define levels for assertion, logging and tracing

#include <scai/lama.hpp>

#include <scai/lama/DenseVector.hpp>
#include <scai/lama/Scalar.hpp>
#include <scai/lama/expression/all.hpp>
#include <scai/lama/matrix/CSRSparseMatrix.hpp>
#include <scai/lama/matutils/MatrixCreator.hpp>

#include <iostream>

using namespace scai::lama;
using namespace std;

int main( int argc, char* argv[] )
{
    File::FileType outType = File::MATRIX_MARKET;

    if ( argc < 2 )
    {
        cout << "Usage: " << argv[0] << " filename[.frm]" << endl;
        return 0;
    }

    string inFileName = argv[1];
    string outFileName = inFileName;

    if ( inFileName.substr( inFileName.size() - 4, 4 ) == ".mtx" )
    {
        // Input file is matrix market, so convert it to binary

        outType = File::BINARY;
        outFileName = inFileName.substr( 0, inFileName.size() - 4 );
        cout << "convert MatrixMarket " << inFileName << " to binary " << outFileName << endl;
    }
    else if ( inFileName.substr( inFileName.size() - 4, 4 ) == ".frm" )
    {
        // Input file is binary, convert it to matrix market
        outFileName = inFileName.substr( 0, inFileName.size() - 4 );
        cout << "convert binary " << inFileName << " to matrix market " << outFileName << endl;
    }

    CSRSparseMatrix<double> m ( inFileName );

    cout << "Read matrix " << inFileName << " : " << m << endl;
    cout << "Write matrix " << outFileName << ", format = " << outType << endl;

    m.writeToFile( outFileName, outType );
}
