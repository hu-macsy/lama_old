/**
 * @file matrix_generator.cpp
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
 * @brief Example program that generates matrices and writes them to a file
 * @author Thomas Brandes
 * @date 15.05.2013
 * @since 1.0.0
 */

// Define levels for assertion, logging and tracing

#include "lama.hpp"

#include <lama/DenseVector.hpp>
#include <lama/Scalar.hpp>
#include <lama/expression/all.hpp>
#include <lama/CommunicatorFactory.hpp>
#include <lama/matrix/CSRSparseMatrix.hpp>
#include <lama/matutils/MatrixCreator.hpp>

#include <iostream>

using namespace lama;
using namespace std;

int main( int argc, char* argv[] )
{
    CommunicatorPtr comm = CommunicatorFactory::get();

    int myRank = comm->getRank();

    const char* filename;
    
    IndexType dimension = 1;
    IndexType stencilType = 3;
    IndexType dimX = 1;
    IndexType dimY = 1;
    IndexType dimZ = 1;

    if ( argc >= 5 )
    {
        filename = argv[1];
        sscanf( argv[2], "%d", &dimension );
        sscanf( argv[3], "%d", &stencilType );
        sscanf( argv[4], "%d", &dimX );
        if ( argc >= 6 ) sscanf( argv[5], "%d", &dimY );
        if ( argc >= 7 ) sscanf( argv[6], "%d", &dimZ );
    }
    else
    {
        if ( myRank == 0 )
        {
            cout << "Usage: " << argv[0] << " <filename> <dim> <stencilType> <dimX> [ <dimY> [ <dimZ> ] ]" << endl;
        }
        return -1;
    }

    cout << "Generate poisson file " << filename << 
            ", dim = " << dimension << ", stencilType = " << stencilType << endl;

    if ( !MatrixCreator<double>::supportedStencilType( dimension, stencilType ) )
    {
        if ( myRank == 0 )
        {
            cout << "Unsupported stencilType " << stencilType << " for dim = " << dimension << endl;
        }
        return -1;
    }

    if ( argc != ( dimension + 4 ) )
    {
        if ( myRank == 0 )
        {
            cout << "Missing values for dim = " << dimension 
                 << ", argc = " << argc << ", expected " << ( dimension + 3 ) << endl;
        }
        return -1;
    }

    cout << "Dimensions: " << dimX;
    if ( dimension > 1 ) cout << " x " << dimY;
    if ( dimension > 2 ) cout << " x " << dimZ;
    cout << endl << endl;

    CSRSparseMatrix<double> m;

    MatrixCreator<double>::buildPoisson( m, dimension, stencilType, dimX, dimY, dimZ );

    DenseVector<double> lhs( m.getDistributionPtr(), 1.0 );
    DenseVector<double> rhs( m * lhs );

    cout << "m = " << m << endl;
    cout << "m has diagonal property = " << m.hasDiagonalProperty() << endl;
    cout << "lhs = " << lhs << endl;
    cout << "rhs = " << rhs << endl;

    m.writeToFile( filename, File::BINARY );
    rhs.writeToFile( filename, File::BINARY );
}

