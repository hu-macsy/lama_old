/**
 * @file lamaGenRandomMatrix.cpp
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
 * @brief Example program that generates matrices and writes them to a file
 * @author Thomas Brandes
 * @date 09.07.2016
 */

// Define levels for assertion, logging and tracing

#include <scai/lama.hpp>

#include <scai/lama/DenseVector.hpp>
#include <scai/lama/Scalar.hpp>
#include <scai/lama/matrix/CSRSparseMatrix.hpp>
#include <scai/lama/matrix/MatrixAssembly.hpp>

#include <scai/common/mepr/TypeListUtils.hpp>
#include <scai/common/Settings.hpp>

#include <iostream>
#include <sstream>

using namespace scai;
using namespace scai::lama;
using namespace scai::dmemo;
using namespace std;

static void replaceInfo( std::string& str, const std::string& info )
{
    size_t i = str.find( "%s" );

    if ( i != string::npos )
    {
        str.replace( i, 2, info );
    }
}

static void printUsage( const char* prog_name )
{
    cout << "Usage: " << prog_name << " <filename> nrows ncols [ fillrate ]" << endl;
    cout << "         filename : name of the output file for matrix, vector" << endl;
    cout << "         fillrate : 1 for dense matrix" << endl;
    cout << "Other command line options:" << endl;
    cout << "     --SCAI_TYPE=float|double|ComplexFloat|ComplexDouble|... specifies value type matrix,vectors" << endl;
    cout << "     --SCAI_IO_BINARY=[0|1]  force formatted or binary output" << endl;
    cout << "Purpose: generate random matrix m, random vector x and b with b = m * x" << endl;
}

/** Get the value type used for this application.
 *
 *  Default value type is DefaultReal( double ) but can be overwritten by
 *
 *   - environment variable SCAI_TYPE=float|double|LongDouble|ComplexFloat|ComplexDouble| ...
 *   - or command line argument --SCAI_TYPE=...
 */

static common::ScalarType getType()
{
    common::ScalarType type = common::TypeTraits<DefaultReal>::stype;

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

template<typename ValueType> 
void generate( const IndexType nrows, const IndexType ncols, const float fillRate, std::string& matrixFileName )
{
    MatrixAssembly<ValueType> assembly;

    for ( IndexType i = 0; i < nrows; ++i )
    {
        for ( IndexType j = 0; j < ncols; ++j )
        {
            bool takeIt = common::Math::randomBool( fillRate );

            if ( takeIt )
            {
                ValueType val = common::Math::random<ValueType>( 1 );
                assembly.push( i, j, val );
            }
        }
    }

    auto m = convert<CSRSparseMatrix<ValueType>>( assembly.buildGlobalCOO( nrows, ncols ) );

    DenseVector<ValueType> x;
    DenseVector<ValueType> b;

    x.setRandom( m.getColDistributionPtr(), 1 );
    b = m * x;

    cout << "m = " << m << endl;
    cout << "m has diagonal property = " << m.hasDiagonalProperty() << endl;
    cout << "x = " << x << endl;
    cout << "b = " << b << endl;
    cout << endl;

    string suffix = FileIO::getSuffix( matrixFileName );

    string vectorXFileName = matrixFileName;
    string vectorBFileName = matrixFileName;

    if ( FileIO::canCreate( suffix ) )
    {
        // known suffix so we can use it directly

        vectorXFileName.replace( vectorXFileName.length() - suffix.length(), 1, "_x." );
        vectorBFileName.replace( vectorBFileName.length() - suffix.length(), 1, "_b." );

        if ( suffix == ".frm" )
        {
            // SAMG format uses two different suffixes for matrix and vector
            // take <filename>.frv instead of <filename>.frm

            vectorXFileName.replace( vectorXFileName.length() - 1, 1, "v" );
            vectorBFileName.replace( vectorBFileName.length() - 1, 1, "v" );
        }
    }
    else
    {
        if ( suffix.length() > 0 )
        {
            cout << "ATTENTION: " << suffix << " is unknown suffix, take SAMG format" << endl;
        }

        matrixFileName += ".frm";
        vectorXFileName += "_x.frv";
        vectorBFileName += "_b.frv";

    }

    cout << "Write matrix to file " << matrixFileName;
    cout << ", vector x to " << vectorXFileName;
    cout << ", and vector b to file " << vectorBFileName << endl;

    b.writeToFile( vectorBFileName );
    cout << "Written vector b (rhs) to file " << vectorBFileName << endl;
    x.writeToFile( vectorXFileName );
    cout << "Written matrix to file " << matrixFileName << endl;
    m.writeToFile( matrixFileName );
    cout << "Written vector x (solution) to file " << vectorXFileName << endl;
}

int main( int argc, const char* argv[] )
{
    common::Settings::parseArgs( argc, argv );

    CommunicatorPtr comm = Communicator::getCommunicatorPtr();

    int myRank = comm->getRank();

    IndexType nrows = 1;
    IndexType ncols = 1;

    std::string matrixFileName;
    float fillRate = 0.3;

    if ( argc >= 4 )
    {
        matrixFileName = argv[1];

        {
            stringstream ss;
            ss << argv[2];
            ss >> nrows;
        }

        {
            stringstream ss;
            ss << argv[3];
            ss >> ncols;
        }

        if ( argc > 4 )
        {
            stringstream ss;
            ss << argv[4];
            ss >> fillRate;
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

    ostringstream info;   // takes nrows x ncols _ fill, fill in %

    info << nrows << "x" << ncols << "_" << static_cast<int>( fillRate * 100 + 0.5f );

    replaceInfo( matrixFileName, info.str() );

    cout << "Generate random file " << matrixFileName << ", fillRate = " << fillRate << endl;

    common::ScalarType stype = getType();

#define DO_GENERATE( ValueType )                                        \
    if ( stype == common::TypeTraits<ValueType>::stype )                \
    {                                                                   \
        generate<ValueType>( nrows, ncols, fillRate, matrixFileName );  \
    }                                                                  

    SCAI_COMMON_LOOP( DO_GENERATE, SCAI_NUMERIC_TYPES_HOST )

#undef DO_GENERATE

    return 0;
}
