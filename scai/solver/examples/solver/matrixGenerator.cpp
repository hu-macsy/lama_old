/**
 * @file matrixGenerator.cpp
 *
 * @license
 * Copyright (c) 2009-2016
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
 * @brief Example program that generates matrices and writes them to a file
 * @author Thomas Brandes
 * @date 15.05.2013
 */

// Define levels for assertion, logging and tracing

#include <scai/lama.hpp>

#include <scai/lama/DenseVector.hpp>
#include <scai/lama/Scalar.hpp>
#include <scai/lama/expression/all.hpp>
#include <scai/lama/StorageIO.hpp>
#include <scai/lama/matrix/CSRSparseMatrix.hpp>
#include <scai/lama/matutils/MatrixCreator.hpp>
#include <scai/common/mepr/TypeListUtils.hpp>

#include <iostream>

using namespace scai;
using namespace scai::lama;
using namespace scai::dmemo;
using namespace std;

void replaceStencil( std::string& str, const std::string& stencil )
{
    size_t i = str.find( "%s" );

    if ( i != string::npos )
    {
        str.replace( i, 2, stencil );
    }
}

void printUsage( const char* prog_name )
{
    cout << "Usage: " << prog_name << " <filename> <dim> <stencilType> <dimX> [ <dimY> [ <dimZ> ] ]" << endl;
    cout << "         filename : name of the output file for matrix, vector" << endl;
    cout << "           filename = <id>.mtx -> generates matrix market format, <id>_v.mtx for vector" << endl;
    cout << "           filename = <id>     -> generates binary format, <id>.frm for matrix, <id>.frv for vector" << endl;
    cout << "           %s in filename is replaced with stencil values, e.g. 2D5P_100_100" << endl;
    cout << "         dim = 1, 2, 3  is dimension of stencil" << endl;
    cout << "         stencilType = 3 (for dim = 1) " << endl;
    cout << "         stencilType = 5, 9 (for dim = 2) " << endl;
    cout << "         stencilType = 7, 19, 27 (for dim = 3) " << endl;
}

/** Define the value type used in this example, take default real type. */

typedef RealType ValueType;

int main( int argc, char* argv[] )
{
    CommunicatorPtr comm = Communicator::getCommunicatorPtr();
    int myRank = comm->getRank();
    std::string matrixFileName;
    IndexType dimension = 1;
    IndexType stencilType = 3;
    IndexType dimX = 1;
    IndexType dimY = 1;
    IndexType dimZ = 1;

    if ( argc >= 5 )
    {
        matrixFileName = argv[1];
        sscanf( argv[2], "%d", &dimension );
        sscanf( argv[3], "%d", &stencilType );
        sscanf( argv[4], "%d", &dimX );

        if ( argc >= 6 )
        {
            sscanf( argv[5], "%d", &dimY );
        }

        if ( argc >= 7 )
        {
            sscanf( argv[6], "%d", &dimZ );
        }
    }
    else
    {
        if ( myRank == 0 )
        {
            printUsage( argv[0] );
        }

        return -1;
    }

    cout << "Generate poisson file " << matrixFileName <<
         ", dim = " << dimension << ", stencilType = " << stencilType << endl;

    if ( !MatrixCreator<ValueType>::supportedStencilType( dimension, stencilType ) )
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

    // Generate name for the stencil
    ostringstream stencilName;
    stencilName << dimension << "D" << stencilType << "P_" << dimX;

    if ( dimension > 1 )
    {
        stencilName << "_" << dimY;
    }

    if ( dimension > 2 )
    {
        stencilName << "_" << dimZ;
    }

    cout << "Stencil is : " << stencilName.str() << endl;
    // replace %s in file name with stencil description
    replaceStencil( matrixFileName, stencilName.str() );
    CSRSparseMatrix<ValueType> m;
    MatrixCreator<ValueType>::buildPoisson( m, dimension, stencilType, dimX, dimY, dimZ );
    DenseVector<ValueType> lhs( m.getRowDistributionPtr(), 1.0 );
    DenseVector<ValueType> rhs( m * lhs );
    cout << "m = " << m << endl;
    cout << "m has diagonal property = " << m.hasDiagonalProperty() << endl;
    cout << "lhs = " << lhs << endl;
    cout << "rhs = " << rhs << endl;
    cout << endl;
    cout << "Solution vector x = ( 1.0, ..., 1.0 ) assumed" << endl;
    cout << "Write matrix and rhs vector to file " << matrixFileName << endl;

    if ( _StorageIO::hasSuffix( matrixFileName, ".mtx" ) )
    {
        std::string vectorFileName = matrixFileName;
        // replace . with _v.
        vectorFileName.replace( vectorFileName.length() - 4, 1, "_v." );
        m.writeToFile( matrixFileName, File::MATRIX_MARKET );
        rhs.writeToFile( vectorFileName, File::MATRIX_MARKET );
        cout << "Written matrix to matrix market file " << matrixFileName  << endl;
        cout << "Written rhs vector to matrix market file " << vectorFileName << endl;
        return 0;
    }
    else 
    {
        std::string vectorFileName = matrixFileName;

        // add suffix frm, frv if not available

        if ( _StorageIO::hasSuffix( matrixFileName, ".frm" ) )
        {
            vectorFileName.replace( vectorFileName.length() - 4, 4, ".frv" );
        }
        else if ( _StorageIO::hasSuffix( matrixFileName, ".frv" ) )
        {
            matrixFileName.replace( matrixFileName.length() - 4, 4, ".frm" );
        }
        else 
        {
            matrixFileName += ".frm";
            vectorFileName += ".frv";
        }

        m.writeToFile( matrixFileName, File::SAMG_FORMAT, common::scalar::INTERNAL, common::scalar::INDEX_TYPE, common::scalar::INDEX_TYPE, true );
        rhs.writeToFile( vectorFileName, File::SAMG_FORMAT, common::scalar::INTERNAL, true );
        cout << "Written matrix to SAMG file " << matrixFileName << endl;
        cout << "Written rhs vector to SAMG file " << vectorFileName << endl;
        return 0;
    }
}
