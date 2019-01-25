/**
 * @file lamaGenStencilMatrix.cpp
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
 * @brief Example program that generates matrices and writes them to a file
 * @author Thomas Brandes
 * @date 15.05.2013
 */

// Define levels for assertion, logging and tracing

#include <scai/lama.hpp>

#include <scai/lama/DenseVector.hpp>
#include <scai/lama/Scalar.hpp>
#include <scai/lama/matrix/CSRSparseMatrix.hpp>
#include <scai/lama/matutils/MatrixCreator.hpp>
#include <scai/lama/io/PartitionIO.hpp>

#include <scai/common/mepr/TypeListUtils.hpp>
#include <scai/common/Settings.hpp>

#include <iostream>
#include <sstream>

using namespace scai;
using namespace scai::lama;
using namespace scai::dmemo;
using namespace std;

#define HOST_PRINT( rank, msg )             \
    {                                           \
        if ( rank == 0 )                        \
        {                                       \
            std::cout << msg << std::endl;      \
        }                                       \
    }                                           \
     
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
    cout << "   other options:" << endl;
    cout << "      --SCAI_IO_TYPE_DATA=float|double|ComplexFloat|ComplexDouble to force output in other format" << endl;
    cout << "      --SCAI_IO_TYPE_INDEX=int|long to fource output of indexes in other format" << endl;
    cout << "      --SCAI_IO_BINARY=0|1 to force formatted or binary format" << endl;
}

struct CommandLineArguments
{
    std::string matrixFileName;

    IndexType dimension;
    IndexType stencilType;
    IndexType dimX;
    IndexType dimY;
    IndexType dimZ;

    CommandLineArguments()
    {
        dimension = 1;
        stencilType = 3;
        dimX = 1;
        dimY = 1;
        dimZ = 1;
    }

    bool parseArguments( PartitionId myRank, int argc, const char* argv[] )
    {
        if ( argc < 5 )
        {
            HOST_PRINT( myRank, "Insufficient arguments." )

            return false;
        }

        matrixFileName = argv[1];
        {
            std::stringstream ss;
            ss << argv[2];
            ss >> dimension;
        }

        {
            std::stringstream ss;
            ss << argv[3];
            ss >> stencilType;
        }

        {
            std::stringstream ss;
            ss << argv[4];
            ss >> dimX;
        }

        if ( argc >= 6 )
        {
            std::stringstream ss;
            ss << argv[5];
            ss >> dimY;
        }

        if ( argc >= 7 )
        {
            std::stringstream ss;
            ss << argv[6];
            ss >> dimZ;
        }

        if ( static_cast<IndexType>( argc ) != ( dimension + 4 ) )
        {
            HOST_PRINT( myRank, "Missing values for dim = " << dimension
                        << ", argc = " << argc << ", expected " << ( dimension + 3 ) );

            return false;
        }

        return true;
    }

    std::string stencilName()
    {
        // Generate name for the stencil

        ostringstream ostream;

        ostream << dimension << "D" << stencilType << "P_" << dimX;

        if ( dimension > 1 )
        {
            ostream << "_" << dimY;
        }

        if ( dimension > 2 )
        {
            ostream << "_" << dimZ;
        }

        return ostream.str();
    }
};

/** Define the value type used in this example, take default real type. */

typedef DefaultReal ValueType;

int main( int argc, const char* argv[] )
{
    common::Settings::parseArgs( argc, argv );

    CommunicatorPtr comm = Communicator::getCommunicatorPtr();
    int myRank = comm->getRank();

    try
    {
        CommandLineArguments cmdArgs;

        bool okay = cmdArgs.parseArguments( myRank, argc, argv );

        if ( !okay )
        {
            if ( myRank == 0 )
            {
                printUsage( argv[0] );
            }

            return -1;
        }

        if ( !MatrixCreator::supportedStencilType( cmdArgs.dimension, cmdArgs.stencilType ) )
        {
            HOST_PRINT( myRank, "Unsupported stencilType " << cmdArgs.stencilType << " for dim = " << cmdArgs.dimension )
            return -1;
        }

        std::string stencilName = cmdArgs.stencilName();

        HOST_PRINT( myRank, "Stencil is : " << stencilName )

        std::string matrixFileName = cmdArgs.matrixFileName;

        // replace %s in file name with stencil description
        replaceStencil( matrixFileName, stencilName );
        CSRSparseMatrix<ValueType> m;

        MatrixCreator::buildPoisson( m, cmdArgs.dimension, cmdArgs.stencilType, cmdArgs.dimX, cmdArgs.dimY, cmdArgs.dimZ );

        auto lhs = denseVectorFill<ValueType>( m.getRowDistributionPtr(), 1 );
        auto rhs = denseVectorEval( m * lhs );

        HOST_PRINT( myRank, "Poisson matrix m = " << m )

        HOST_PRINT( myRank, "lhs = " << lhs << ", is all 1.0" )
        HOST_PRINT( myRank, "rhs = " << rhs << ", is m * lhs" )

        std::string suffix = FileIO::getSuffix( matrixFileName );

        std::string vectorFileName = matrixFileName;

        if ( FileIO::canCreate( suffix ) )
        {
            // known suffix so we can use it directly

            if ( suffix == ".frm" )
            {
                // SAMG format uses two different suffixes for matrix and vector
                // take <filename>.frv instead of <filename>.frm

                vectorFileName.replace( vectorFileName.length() - 4, 4, ".frv" );
            }
            else
            {
                // take <filename>_v.<suffix> for <filename>.<suffix>

                vectorFileName.replace( vectorFileName.length() - suffix.length(), 1, "_v." );
            }
        }
        else
        {
            if ( suffix.length() > 0 )
            {
                HOST_PRINT( myRank, "ATTENTION: " << suffix << " is unknown suffix, take SAMG format" )
            }

            matrixFileName += ".frm";
            vectorFileName += ".frv";

        }

        HOST_PRINT( myRank, "Write matrix to file " << matrixFileName << " and rhs vector to file " << vectorFileName )

        if ( PartitionIO::isPartitionFileName( matrixFileName ) )
        {
            // check if is necessary to write the mapping

            IndexType nb = m.getRowDistribution().getBlockDistributionSize();

            if ( nb == invalidIndex )
            {
                HOST_PRINT( myRank, "WARNING: matrix has no block distribution" )
            }

            std::string distFileName;

            if ( common::Settings::getEnvironment( distFileName, "SCAI_DISTRIBUTION" ) )
            {
                PartitionIO::write( m.getRowDistribution(), distFileName );

                HOST_PRINT( myRank, "written distribution to " << distFileName )
            }
        }

        m.writeToFile( matrixFileName );
        rhs.writeToFile( vectorFileName );

        HOST_PRINT( myRank, "Written matrix to file " << matrixFileName )
        HOST_PRINT( myRank, "Written rhs vector to file " << vectorFileName )

    }
    catch ( common::Exception& e )
    {
        HOST_PRINT( myRank, "Caught exception: " << e.what() << "\n\nTERMINATE due to error" )
        return -1;
    }

    return 0;
}
