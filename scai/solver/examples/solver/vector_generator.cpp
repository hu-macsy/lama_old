/**
 * @file vector_generator.cpp
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
 * @brief Program to generate files containing vectors
 * @author Thomas Brandes
 * @date 14.01.2016
 */

#include <scai/lama.hpp>

#include <scai/lama/DenseVector.hpp>
#include <scai/lama/Scalar.hpp>
#include <scai/lama/matrix/Matrix.hpp>
#include <scai/lama/StorageIO.hpp>
#include <scai/dmemo/NoDistribution.hpp>
#include <scai/lama/expression/all.hpp>

#include <iostream>
#include <algorithm>

using namespace scai;
using namespace lama;
using namespace dmemo;
using namespace std;

static bool isNumber( const char* arg )
{
    int len = strlen( arg );

    for ( int i = 0; i < len; ++i )
    {
        if ( isdigit( arg[i] ) )
        {
            continue;
        }
        return false;
    }
    return true;
}

static bool isValue( const char* arg )
{
    int len = strlen( arg );

    for ( int i = 0; i < len; ++i )
    {
        if ( isdigit( arg[i] ) )
        {
            continue;
        }
        if ( arg[i] == '.' || arg[i] == ',' || arg[i] == '-' || arg[i] == ' ' )
        {
            continue;
        }
        return false;
    }
    return true;
}

struct CommandLineOptions
{
    string outFileName;
    string matFileName;

    File::FileType outFileType;
    common::scalar::ScalarType outDataType;

    Scalar value;   // value for the vector

    bool random;    // if true generate random numbers

    IndexType size;

    CommandLineOptions()
    {
        outFileName = "";
        matFileName = "";
        outFileType = File::DEFAULT;
        outDataType = common::scalar::INTERNAL;       // same as input data type
        value       = Scalar( 1 );
        size        = 0;
        random      = false;
    }

    bool parseOption( const std::string& option )
    {
        if ( option == "-a" )
        {
            outFileType = File::FORMATTED;
        }
        else if ( option == "-mm" )
        {
            outFileType = File::MATRIX_MARKET;
        }
        else if ( option == "-b" )
        {
            outFileType = File::BINARY;
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
        else if ( option == "-random" )
        {
            random = true;
        } 
        else if ( isNumber( option.c_str() ) && size == 0 )
        {
            std::istringstream is( option );
            is >> size;
        }
        else if ( isValue( option.c_str() ) )
        {
            ComplexDouble x;
            std::string vstr = option;
            std::replace( vstr.begin(), vstr.end(), ',', ' '); // replace all ',' to ' '
            std::istringstream is( vstr );
            is >> x;
            cout << "Read value from " << vstr << ", value = " << x << endl;
            value = x;
        }
        else if ( outFileName == "" )
        {    
            outFileName = option;
        }
        else if ( matFileName == "" && _StorageIO::fileExists( option ) )
        {
            matFileName = option;
        }
        else
        {
            return false;  // not recognized
        }
        return true;
    }

    void checkOutFileName()
    {
        if ( outFileName == "" )
        {
            cout << "No outfile name specified, take 'vector'" << endl;
            outFileName = "vector";
        }
    }

    void checkOutFileType()
    {
        if ( outFileType != File::DEFAULT ) return;

        if ( _StorageIO::hasSuffix( outFileName, ".mtx" ) )
        {
             outFileType = File::MATRIX_MARKET;
        }
        else
        {
             outFileType = File::BINARY;
        }

        cout << "No output file type specified, take " << outFileType << endl;
    }

    void checkOutDataType()
    {
        if ( outDataType != common::scalar::INTERNAL ) return;

        // take double or complex double

        if ( conj( value ) == value )
        {
            outDataType = common::scalar::DOUBLE;
        }
        else
        {
            outDataType = common::scalar::DOUBLE_COMPLEX;
        }
       
        cout << "No output data type specified, take " << outDataType << " due to value = " << value << endl;
    }
};

void printUsage( const char* progName )
{
    cout << "Usage: " << progName << " [-b|-a|-mm] [-s|-d|-c|-z] outfile_name <size> <val> [matrix_filename]" << endl;
    cout << "   outfile_name is filename for vector output" << endl;
    cout << "    -b generates binary file [default if filename has not suffix .mtx]" << endl;
    cout << "    -a generates formatted file" << endl;
    cout << "    -mm generates matrix market format[ default if filename has suffix .mtx]" << endl;
    cout << "   size is the number of elements in the vector" << endl;
    cout << "   val is the value for each entry" << endl;
    cout << "    -random each entry is multiplied with a random value from 0..1" << endl;
    cout << "    -d output format is double precision [default]" << endl;
    cout << "    -c output format is complex" << endl;
    cout << "    -z output format is double complex" << endl;
    cout << "   matrix_filename : if set, compute vector as rhs of matrix * vector" << endl;
}

int main( int argc, char* argv[] )
{
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
            cout << endl;
            cout << "Option " << argv[i] << " not recognized" << endl;
            cout << endl;
            printUsage( argv[0] );
            return -1;
        }
    }

    options.checkOutFileName();
    options.checkOutFileType();
    options.checkOutDataType();

    cout << "Generate vector ( size = " << options.size << ", val = " << options.value << " )" << endl;

    // use vector of outDataType so no information is lost

    common::shared_ptr<Matrix> matrix;

    VectorCreateKeyType vectorType( DENSE, options.outDataType );
    common::shared_ptr<Vector> v ( Vector::create( vectorType ) );

    if ( options.matFileName != "" )
    { 
        MatrixCreateKeyType matrixType( Format::CSR, options.outDataType );
        matrix.reset( Matrix::create( MatrixCreateKeyType ( matrixType ) ) );
        matrix->readFromFile( options.matFileName );
        cout << "Read in matrix from file " << options.matFileName << ": " << *matrix << endl;
    }

    if ( options.size == 0 && matrix.get() )
    {
        v->resize( matrix->getColDistributionPtr() );
    }
    else
    {
        DistributionPtr dist ( new NoDistribution( options.size ) );
        v->resize( dist );
    }

    cout << "Vector (uninitialized): " << *v << endl;

    *v = options.value;

    if ( options.random )
    {
        using namespace hmemo;

        // we know what we do here, so const_cast is okay

        _HArray& vLocal = const_cast<_HArray&>( v->getLocalValues() );

        ContextPtr host = Context::getHostPtr();

        HArray<double> randomValues;

        {
            IndexType n = vLocal.size();

            std::srand( 171451 );

            WriteOnlyAccess<double> write( randomValues, host, n );

            for ( int i = 0; i < n; ++i )
            {
                write[i] = static_cast<double>( rand() ) / static_cast<double>( RAND_MAX );
            }
        }

        // multiply random numbers with the scalar value

        HArrayUtils::set( vLocal, randomValues, common::reduction::MULT, host );
    }

    cout << "Vector generated: " << *v << endl;

    if ( matrix.get() )
    { 
        VectorCreateKeyType vectorType( DENSE, options.outDataType );
        common::shared_ptr<Vector> rhs ( Vector::create( vectorType ) );
        *rhs = *matrix * *v;
        v = rhs;
        cout << "Vector now rhs of multiplication with matrix: " << *v << endl;
    }

    cout << "write to output file " << options.outFileName;
    cout << ", format = " << options.outFileType;
    cout << ", data type = " << options.outDataType;
    cout << endl;

    v->writeToFile( options.outFileName, options.outFileType, common::scalar::INTERNAL );

    cout << "Done." << endl;
}
