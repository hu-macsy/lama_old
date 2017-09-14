/**
 * @file PoissonInputSetCreator.cpp
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
 * @brief PoissonInputSetCreator.cpp
 * @author Jiri Kraus
 * @date 06.05.2010
 */

#include <scai/lama/benchmark/PoissonInputSetCreator.hpp>

#include <scai/lama/matrix/CSRSparseMatrix.hpp>

#include <scai/dmemo/BlockDistribution.hpp>
#include <scai/dmemo/GeneralDistribution.hpp>
#include <scai/lama/matutils/MatrixCreator.hpp>
#include <scai/lama/expression/MatrixVectorExpressions.hpp>

#include <sstream>
#include <iostream>
#include <algorithm>
#include <functional>
#include <vector>
#include <ctime>
#include <cstdlib>
#include <cstdio>

namespace
{

std::vector<std::string> split( const std::string& params, const char seperator )
{
    std::vector<std::string> args;
    size_t found = std::string::npos;

    do
    {
        size_t prevFound = found + 1;
        found = params.find( seperator, prevFound );
        args.push_back( params.substr( prevFound, found - prevFound ) );
    }
    while ( found != std::string::npos );

    return args;
}
}

using namespace scai;

SCAI_LOG_DEF_LOGGER( PoissonInputSetCreator::logger, "InputSetCreator.PoissonInputSetCreator" );

PoissonInputSetCreator::PoissonInputSetCreator()
{
    SCAI_LOG_INFO( logger, "PoissonInputSetCreator" );
}

PoissonInputSetCreator::~PoissonInputSetCreator()
{
    SCAI_LOG_INFO( logger, "~PoissonInputSetCreator" );
}

PoissonInputSetCreator::InputSetType* PoissonInputSetCreator::create() const
{
    std::ostringstream message;
    message << "PoissonInputSetCreator::create: This InputSetCreator requires the following arguments:" << std::endl
            << "    aDbP_dimX_dimY_dimZ," << std::endl
            << "where a is in (1, 2, 3) and b is in (3, 5, 7, 9, 19, 27) and dimY "
            << "and dimZ are needed in case of unequality to one, only." << std::endl;
    throw bf::BFException( message.str() );
}

PoissonInputSetCreator::InputSetType* PoissonInputSetCreator::createSet( const std::string& params )
{
    SCAI_LOG_INFO( logger, "create, params = " << params );

    // 2D5P_10000_10
    std::vector<std::string> args = split( params, '_' );

    if ( args.size() < 2 || args.size() > 4 )
    {
        std::ostringstream message;
        message << "PoissonInputSetCreator::create( ): Wrong number of " << "parameters (" << args.size()
                << "). Expected:" << std::endl << "    aDbP_dimX_dimY_dimZ," << std::endl
                << "where a is in (1, 2, 3) and b is in (3, 5, 7, 9, 19, 27) and dimY "
                << "and dimZ are needed in case of unequality to one, only." << std::endl;
        throw bf::BFException( message.str() );
    }

    IndexType stencilType;
    IndexType dimension;
    IndexType dimX;
    IndexType dimY = 1;
    IndexType dimZ = 1;
    std::vector<std::string> vals;

    if ( args[0].find( 'D', 0 ) == std::string::npos || args[0].find( 'P', 0 ) == std::string::npos )
    {
        std::ostringstream message;
        message << "PoissonInputSetCreator::create( ): Expected 'D' and 'P' within " << args[0]
                << " but at least one was not found!" << std::endl;
        throw bf::BFException( message.str() );
    }

    vals = split( args[0], 'D' );

    if ( vals.size() != 2 )
    {
        std::ostringstream message;
        message << "PoissonInputSetCreator::create( ): Illegal format for " << args[0] << std::endl
                << "Should have been aDbP, with a=1,2,3 and b=3,5,7,9,19 or 27." << std::endl;
        throw bf::BFException( message.str() );
    }

    vals[1] = split( vals[1], 'P' )[0];
    std::istringstream convert( vals[0] );
    convert >> dimension;

    SCAI_LOG_DEBUG( logger, "Stencil dimension = " << dimension );

    if ( convert.fail() )
    {
        std::ostringstream message;
        message << "PoissonInputSetCreator::create( ): Unable to convert '" << convert.str()
                << "' (should have been dimension)." << std::endl;
        throw bf::BFException( message.str() );
    }

    convert.clear();
    convert.str( vals[1] );
    convert >> stencilType;

    if ( convert.fail() )
    {
        std::ostringstream message;
        message << "PoissonInputSetCreator::create( ): Unable to convert '" << convert.str()
                << "' (should have been stencilType)." << std::endl;
        throw bf::BFException( message.str() );
    }

    convert.clear();

    if ( dimension != ( IndexType ) args.size() - 1 )
    {
        std::ostringstream message;
        message << "PoissonInputSetCreator::create( ): '" << params << "' said " << "to have " << dimension
                << " dimensions, but had " << ( args.size() - 1 ) << " dimension parameters." << std::endl;
        throw bf::BFException( message.str() );
    }

    switch ( dimension )
    {
        case 3:
            convert.str( args[3] );
            convert >> std::dec >> dimZ;

            if ( convert.fail() )
            {
                std::ostringstream message;
                message << "PoissonInputSetCreator::create( ): Unable to convert '" << convert.str()
                        << "' (should have been dimZ)." << std::endl;
                throw bf::BFException( message.str() );
            }

            convert.clear();
        case 2:
            convert.str( args[2] );
            convert >> std::dec >> dimY;

            if ( convert.fail() )
            {
                std::ostringstream message;
                message << "PoissonInputSetCreator::create( ): Unable to convert '" << convert.str()
                        << "' (should have been dimY)." << std::endl;
                throw bf::BFException( message.str() );
            }

            convert.clear();
        case 1:
            convert.str( args[1] );
            convert >> std::dec >> dimX;

            if ( convert.fail() )
            {
                std::ostringstream message;
                message << "PoissonInputSetCreator::create( ): Unable to convert '" << convert.str()
                        << "' (should have been dimX)." << std::endl;
                throw bf::BFException( message.str() );
            }

            convert.clear();
            break;
        default:
            std::ostringstream message;
            message << "PoissonInputSetCreator::create( ): Invalid number of " << "dimensions (" << dimension
                    << "). Expected 1, 2 or 3." << std::endl;
            throw bf::BFException( message.str() );
    }

    SCAI_LOG_INFO( logger, "dimX = " << dimX << ", dimY = " << dimY << ", dimZ = " << dimZ );

    bool exception = false;

    if ( dimZ == 1 && dimY == 1 )
    {
        exception = ( stencilType != 3 );
    }
    else if ( dimZ == 1 )
    {
        exception = ( stencilType != 5 && stencilType != 9 );
    }
    else if ( dimZ >= 1 && dimY >= 1 )
    {
        exception = ( stencilType != 7 && stencilType != 19 && stencilType != 27 );
    }
    else
    {
        exception = true;
    }

    if ( exception )
    {
        std::ostringstream message;
        message << "PoissonInputSetCreator::create( ): Wrong values for parameters:" << std::endl << "Valid values are:"
                << std::endl << " dimY  | dimZ  | stencilType " << std::endl << "-------+-------+-------------"
                << std::endl << "  == 1 |  == 1 |           3 " << std::endl << "  >  1 |  == 1 |        5, 9 "
                << std::endl << "  >= 1 |  >= 1 |   7, 19, 27 " << std::endl;
        throw bf::BFException( message.str() );
    }

    common::unique_ptr<lama::CSRSparseMatrix<double> > matrixA( new lama::CSRSparseMatrix<double>() );

    lama::MatrixCreator::buildPoisson( *matrixA, dimension, stencilType, dimX, dimY, dimZ );

    dmemo::DistributionPtr rowDistribution = matrixA->getRowDistributionPtr();
    dmemo::DistributionPtr colDistribution = matrixA->getColDistributionPtr();

    common::unique_ptr<lama::DenseVector<double> > vectorX( new lama::DenseVector<double>( colDistribution, 1.0 ) );
    common::unique_ptr<lama::DenseVector<double> > vectorY( new lama::DenseVector<double>( rowDistribution, 0.0 ) );

    SCAI_LOG_DEBUG( logger, "distributed vectors X, Y allocated: " << *vectorX << ", " << *vectorY );

    // Note: in previous versions Y was set when building Poisson matrices, now we do just mult

    *vectorY = ( *matrixA ) * ( *vectorX );

    SCAI_LOG_INFO( logger, "created input set: A = " << *matrixA << ", X = " << *vectorX << ", Y = " << *vectorY );

    return new InputSetType( id(), 1.0, 1.0, vectorX, vectorY, matrixA );
}

PoissonInputSetCreator::InputSetType* PoissonInputSetCreator::create( const std::string& params ) const
{
    return createSet( params );
}

const std::string& PoissonInputSetCreator::id()
{
    static const std::string id = "Poisson";
    return id;
}

const std::string& PoissonInputSetCreator::getId() const
{
    return id();
}

LAMA_INPUTSET_REGISTRATION( PoissonInputSetCreator );
