/**
 * @file lamaGenVector.cpp
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
 * @brief Program to generate files containing vectors
 * @author Thomas Brandes
 * @date 14.01.2016
 */

#include <scai/lama.hpp>

#include <scai/lama/DenseVector.hpp>
#include <scai/lama/matrix/_Matrix.hpp>
#include <scai/dmemo/NoDistribution.hpp>

#include <scai/common/Settings.hpp>

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

static common::ScalarType getType()
{
    common::ScalarType type = common::TypeTraits<double>::stype;

    std::string val;

    if ( scai::common::Settings::getEnvironment( val, "SCAI_TYPE" ) )
    {
        scai::common::ScalarType env_type = scai::common::str2ScalarType( val.c_str() );

        if ( env_type == scai::common::ScalarType::UNKNOWN )
        {
            std::cout << "SCAI_TYPE=" << val << " illegal, is not a scalar type" << std::endl;
        }

        type = env_type;
    }

    return type;
}

struct CommandLineOptions
{
    string outFileName;
    string matFileName;

    common::ScalarType outDataType;

    intern::Scalar value;   // value for the vector

    bool random;    // if true generate random numbers

    IndexType size;

    CommandLineOptions()
    {
        outFileName = "";
        matFileName = "";
        outDataType = getType();
        value       = 1;
        size        = 0;
        random      = false;
    }

    bool parseOption( const std::string& option )
    {
        if ( option == "-random" )
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
#ifdef SCAI_COMPLEX_SUPPORTED
            ComplexDouble x;
#else
            double x;
#endif
            std::string vstr = option;
            std::replace( vstr.begin(), vstr.end(), ',', ' ' ); // replace all ',' to ' '
            std::istringstream is( vstr );
            is >> x;
            cout << "Read value from " << vstr << ", value = " << x << endl;
            value = x;
        }
        else if ( outFileName == "" )
        {
            outFileName = option;
        }
        else if ( matFileName == "" && FileIO::fileExists( option ) )
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

    void checkOutDataType()
    {
        if ( outDataType != common::ScalarType::INTERNAL )
        {
            return;
        }

        // take double or complex double

        outDataType = common::ScalarType::DOUBLE;

        cout << "No output data type specified, take " << outDataType << endl;
    }
};

void printUsage( const char* progName )
{
    cout << "Usage: " << progName << " [--SCAI_var=val] outfile_name <size> <val> [matrix_filename]" << endl;
    cout << "   outfile_name is filename for vector output" << endl;
    cout << "    --SCAI_IO_BINARY=0|1 force formatted or binary output" << endl;
    cout << "    --SCAI_TYPE=float|double|LongDouble|ComplexFloat|ComplexDouble|ComplexLong value type" << endl;
    cout << "   size is the number of elements in the vector" << endl;
    cout << "   val is the value for each entry" << endl;
    cout << "    -random each entry is multiplied with a random value from 0..1" << endl;
    cout << "   matrix_filename : if set, compute vector as rhs of matrix * vector" << endl;
}

template<typename ValueType>
void generate( const CommandLineOptions& options )
{
    // use vector of outDataType so no information is lost

    CSRSparseMatrix<ValueType> matrix;
    DenseVector<ValueType> v;
 
    ValueType initValue = options.value.getValue<ValueType>();

    if ( options.matFileName != "" )
    {
        matrix.readFromFile( options.matFileName );
        cout << "Read in matrix from file " << options.matFileName << ": " << matrix << endl;

        v.setSameValue( matrix.getColDistributionPtr(), initValue );
    }
    else
    {
        v.setSameValue( options.size, initValue );
    }

    cout << "Vector (initialized): " << v << endl;

    if ( options.random )
    {
        // generate random number for the vector, in range (0, 1)

        v.fillRandom( 1 ); 

        // scale random numbers from 0 .. 1 with options.value

        v *= initValue;
    }

    cout << "Vector generated: " << v << endl;

    if ( v.getDistribution() == matrix.getColDistribution() )
    {
        DenseVector<ValueType> rhs;
        rhs = matrix * v;
        v = rhs;
        cout << "Vector now rhs of multiplication with matrix: " << v << endl;
    }

    cout << "write to output file " << options.outFileName << ", data type = " << options.outDataType;
    cout << endl;
    v.writeToFile( options.outFileName );

    cout << "Done." << endl;
}

int main( int argc, const char* argv[] )
{
    common::Settings::parseArgs( argc, argv );

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
    options.checkOutDataType();
    cout << "Generate vector ( size = " << options.size << ", val = " << options.value << " )" << endl;

#define DO_GENERATE( ValueType )                                        \
    if ( options.outDataType == common::TypeTraits<ValueType>::stype )  \
    {                                                                   \
        generate<ValueType>( options );                                 \
    }                                                                  

    SCAI_COMMON_LOOP( DO_GENERATE, SCAI_NUMERIC_TYPES_HOST )

}
