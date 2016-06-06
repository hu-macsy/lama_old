/**
 * @file solver/examples/solver/matrix_convert.cpp
 *
 * @license
 * Copyright (c) 2009-2016
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the Library of Accelerated Math Applications (LAMA).
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
#include <scai/lama/StorageIO.hpp>
#include <scai/lama/expression/all.hpp>
#include <scai/lama/matrix/CSRSparseMatrix.hpp>
#include <scai/lama/matutils/MatrixCreator.hpp>

#include <iostream>

using namespace scai;
using namespace scai::lama;
using namespace std;

struct CommandLineOptions
{
    string inFileName;
    string outFileName;

    File::FileType outFileType;
    bool writeBinary;
    common::scalar::ScalarType inDataType;
    common::scalar::ScalarType outDataType;

    CommandLineOptions()
    {
        inFileName = "";
        outFileName = "";
        writeBinary = false;
        outFileType = File::DEFAULT;
        inDataType  = common::scalar::UNKNOWN;        // needs to be determined
        outDataType = common::scalar::INTERNAL;       // same as input data type
    }

    bool parseOption( const std::string& option )
    {
        if ( option == "-a" )
        {
            outFileType = File::SAMG_FORMAT;
            writeBinary = false;
        }
        else if ( option == "-mm" )
        {
            outFileType = File::MATRIX_MARKET;
            writeBinary = false;
        }
        else if ( option == "-b" )
        {
            outFileType = File::SAMG_FORMAT;
            writeBinary = true;
        }
        else if ( option == "-S" )
        {
            inDataType = common::scalar::FLOAT;
        }
        else if ( option == "-C" )
        {
            inDataType = common::scalar::COMPLEX;
        }
        else if ( option == "-D" )
        {
            inDataType = common::scalar::DOUBLE;
        }
        else if ( option == "-Z" )
        {
            inDataType = common::scalar::DOUBLE_COMPLEX;
        }
        else if ( option == "-s" )
        {
            outDataType = common::scalar::FLOAT;
        }
        else if ( option == "-c" )
        {
            outDataType = common::scalar::COMPLEX;
        }
        else if ( option == "-d" )
        {
            outDataType = common::scalar::DOUBLE;
        }
        else if ( option == "-z" )
        {
            outDataType = common::scalar::DOUBLE_COMPLEX;
        }
        else if ( inFileName == "" )
        {
            inFileName = option;
        }
        else if ( outFileName == "" )
        {
            outFileName = option;
        }
        else
        {
            return false;  // not recognized
        }

        return true;
    }

    void checkInFileName()
    {
        if ( inFileName == "" )
        {
            COMMON_THROWEXCEPTION( "No input filename specified" )
        }
    }

    void checkOutFileName()
    {
        if ( outFileName != "" ) return;

        size_t i = inFileName.find( "." );

        if ( i == string::npos )
        {
            outFileName = inFileName;
        }
        else
        {
            outFileName = inFileName.substr( 0, i );
        }

        cout << "outFileName set to " << outFileName << endl;
    }

    /** @brief determine outfile type if not specified */

    void checkOutFileType()
    {
        if ( outFileType != File::DEFAULT ) return;

        if ( _StorageIO::hasSuffix( inFileName, ".mtx" ) )
        {
             outFileType = File::SAMG_FORMAT;
        }
        else
        {
            outFileType = File::MATRIX_MARKET;
        }
    }
};

void convertMatrix(
    const std::string& inFileName,
    const common::scalar::ScalarType inDataType,
    const std::string& outFileName,
    const File::FileType outFileType,
    const common::scalar::ScalarType outDataType,
    const bool writeBinary )
{
    MatrixCreateKeyType matrixType( Format::CSR, inDataType );
    common::shared_ptr<Matrix> m ( Matrix::create( matrixType ) );

    m->readFromFile( inFileName );

    cout << "read matrix from " << inFileName << " : " << *m << endl;
    cout << "write matrix to " << outFileName << ", format = " << outFileType << ", type = " << outDataType << endl;

    m->writeToFile( outFileName, outFileType, outDataType, common::scalar::INDEX_TYPE, common::scalar::INDEX_TYPE, writeBinary );
}

void convertVector(
    const std::string& inFileName,
    const common::scalar::ScalarType inDataType,
    const std::string& outFileName,
    const File::FileType outFileType,
    const common::scalar::ScalarType outDataType,
    const bool writeBinary )
{
    // Note: inFileType is given implicitly by the input file
    // use vector of inDataType so no information is lost

    VectorCreateKeyType vectorType( Vector::DENSE, inDataType );
    common::shared_ptr<Vector> v ( Vector::create( vectorType ) );

    v->readFromFile( inFileName );

    cout << "read vector from " << inFileName << " : " << *v << endl;
    cout << "write vector to " << outFileName << ", format = " << outFileType << ", type = " << outDataType << endl;

    v->writeToFile( outFileName, outFileType, outDataType, writeBinary );
}

void printUsage( const char* progName )
{
    cout << "Usage: " << progName << " [-b|-a|-mm] infile_name [outfile_name]" << endl;
    cout << "  infile_name : name of file with input matrix or vector" << endl;
    cout << "    -S input format is single precision" << endl;
    cout << "    -D input format is double precision" << endl;
    cout << "    -C input format is complex" << endl;
    cout << "    -Z input format is double complex" << endl;
    cout << "    -s output format is single precision" << endl;
    cout << "    -d output format is double precision" << endl;
    cout << "    -c output format is complex" << endl;
    cout << "    -z output format is double complex" << endl;
    cout << "    -b converts to binary SAMG file" << endl;
    cout << "    -a converts to ascii SAMG file" << endl;
    cout << "    -mm converts to matrix market" << endl;
}

int main( int argc, char* argv[] )
{
    bool isVector = false;

    CommandLineOptions options;

    if ( argc < 2 )
    {
        printUsage( argv[0] );
        return 0;
    }

    for ( int i = 1; i < argc; ++i )
    {
        bool done = options.parseOption( argv[i] );

        if ( !done )
        {
            printUsage( argv[0] );
        }
    }

    options.checkInFileName();
    options.checkOutFileName();

    options.checkOutFileType();

    if ( _StorageIO::hasSuffix( options.inFileName, ".mtx" ) )
    {
        bool isSym, isPat;
        IndexType m, n, nz;

        lama::FileStream inFile( options.inFileName, std::ios::in );

        _StorageIO::readMMHeader( m, n, nz, isSym, isPat, inFile );
        cout << "MM header of " << options.inFileName << ": m = " << m << ", n = " << n << ", nz = " << nz << endl;
        isVector = n == 1;
    }
    else if ( _StorageIO::hasSuffix( options.inFileName, ".frv" ) )
    {
        isVector = true;
    }

    cout << "convert " << ( isVector ? "vector" : "matrix" );
    cout << " " << options.inFileName << " -> " << options.outFileName ;
    cout << ", outFileType = " << options.outFileType;
    cout << ", inDataType = " << options.inDataType;
    cout << ", outDataType = " << options.outDataType << endl;

    if ( options.inDataType == common::scalar::UNKNOWN )
    {
        cout << "No input data type specified, take Double as default" << endl;
        options.inDataType = common::scalar::DOUBLE;
    }

    if ( isVector )
    {
        convertVector( options.inFileName, options.inDataType, options.outFileName, options.outFileType, options.outDataType, options.writeBinary );
    }
    else
    {
        convertMatrix( options.inFileName, options.inDataType, options.outFileName, options.outFileType, options.outDataType, options.writeBinary );
    }
}
