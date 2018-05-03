/**
 * @file PoissonInputSet.cpp
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
 * @brief PoissonInputSet.cpp
 * @author Thomas Brandes
 * @date 15.09.2017
 */

#include <scai/lama/benchmark/PoissonInputSet.hpp>
#include <scai/lama/expression/MatrixVectorExpressions.hpp>
#include <scai/lama/matutils/MatrixCreator.hpp>

#include <scai/common/Settings.hpp>

namespace scai
{

namespace lama
{

void PoissonInputSet::parsePoissonSpecification()
{
    std::vector<std::string> args;   // tokens of 2D5P_10000_10

    common::Settings::tokenize( args, mArgument, "_" );

    if ( args.size() < 2 || args.size() > 4 )
    {
        std::ostringstream message;
        message << "PoissonInputSetCreator::create( ): Wrong number of " << "parameters (" << args.size()
                << "). Expected:" << std::endl << "    aDbP_sizeX_sizeY_sizeZ," << std::endl
                << "where a is in (1, 2, 3) and b is in (3, 5, 7, 9, 19, 27) and sizeY "
                << "and dimZ are needed in case of unequality to one, only." << std::endl;
        throw benchmark::BFException( message.str() );
    }

    std::vector<std::string> vals;

    if ( args[0].find( 'D', 0 ) == std::string::npos || args[0].find( 'P', 0 ) == std::string::npos )
    {
        std::ostringstream message;
        message << "PoissonInputSetCreator::create( ): Expected 'D' and 'P' within " << args[0]
                << " but at least one was not found!" << std::endl;
        throw benchmark::BFException( message.str() );
    }

    common::Settings::tokenize( vals, args[0], "D" );

    if ( vals.size() != 2 )
    {
        std::ostringstream message;
        message << "PoissonInputSetCreator::create( ): Illegal format for " << args[0] << std::endl
                << "Should have been aDbP, with a=1,2,3 and b=3,5,7,9,19 or 27." << std::endl;
        throw benchmark::BFException( message.str() );
    }

    std::istringstream convert( vals[0] );
    convert >> mDimension;

    SCAI_LOG_DEBUG( logger, "Stencil dimension = " << mDimension );

    if ( convert.fail() )
    {
        std::ostringstream message;
        message << "PoissonInputSetCreator::create( ): Unable to convert '" << convert.str()
                << "' (should have been dimension)." << std::endl;
        throw benchmark::BFException( message.str() );
    }

    common::Settings::tokenize( vals, vals[1], "P" );

    convert.clear();
    convert.str( vals[0] );
    convert >> mStencilType;

    if ( convert.fail() )
    {
        std::ostringstream message;
        message << "PoissonInputSetCreator::create( ): Unable to convert '" << convert.str()
                << "' (should have been stencilType)." << std::endl;
        throw benchmark::BFException( message.str() );
    }

    convert.clear();

    if ( mDimension != ( IndexType ) args.size() - 1 )
    {
        std::ostringstream message;
        message << "PoissonInputSetCreator::create( ): '" << mArgument << "' said " << "to have " << mDimension
                << " dimensions, but had " << ( args.size() - 1 ) << " dimension parameters." << std::endl;
        throw benchmark::BFException( message.str() );
    }

    switch ( mDimension )
    {
        case 3:
            convert.str( args[3] );
            convert >> std::dec >> mSizeZ;

            if ( convert.fail() )
            {
                std::ostringstream message;
                message << "PoissonInputSetCreator::create( ): Unable to convert '" << convert.str()
                        << "' (should have been dimZ)." << std::endl;
                throw benchmark::BFException( message.str() );
            }

            convert.clear();
        case 2:
            convert.str( args[2] );
            convert >> std::dec >> mSizeY;

            if ( convert.fail() )
            {
                std::ostringstream message;
                message << "PoissonInputSetCreator::create( ): Unable to convert '" << convert.str()
                        << "' (should have been dimY)." << std::endl;
                throw benchmark::BFException( message.str() );
            }

            convert.clear();
        case 1:
            convert.str( args[1] );
            convert >> std::dec >> mSizeX;

            if ( convert.fail() )
            {
                std::ostringstream message;
                message << "PoissonInputSetCreator::create( ): Unable to convert '" << convert.str()
                        << "' (should have been dimX)." << std::endl;
                throw benchmark::BFException( message.str() );
            }

            convert.clear();
            break;
        default:
            std::ostringstream message;
            message << "PoissonInputSetCreator::create( ): Invalid number of " << "dimensions (" << mDimension
                    << "). Expected 1, 2 or 3." << std::endl;
            throw benchmark::BFException( message.str() );
    }

}

PoissonInputSet::PoissonInputSet( const std::string argument ) : 

    LAMAInputSet(),
    mArgument( argument ),
    mDimension( 1 ),
    mStencilType( 3 ),
    mSizeX( 100 ),
    mSizeY(  1  ),
    mSizeZ(  1  )

{
    SCAI_LOG_INFO( logger, "create, argument = " << argument );

    if ( argument == "" )
    {
        mArgument = "2D5P_100_150";
    }

    parsePoissonSpecification();

    SCAI_LOG_INFO( logger, "Stencil: #dim = " << mDimension << ", type = " << mStencilType
                            << ", size = " << mSizeX << " x " << mSizeY << " x " << mSizeZ );

    bool exception = false;

    if ( mSizeZ == 1 && mSizeY == 1 )
    {
        exception = ( mStencilType != 3 );
    }
    else if ( mSizeZ == 1 )
    {
        exception = ( mStencilType != 5 && mStencilType != 9 );
    }
    else if ( mSizeZ >= 1 && mSizeY >= 1 )
    {
        exception = ( mStencilType != 7 && mStencilType != 19 && mStencilType != 27 );
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
        throw benchmark::BFException( message.str() );
    }

    mAlpha = 1.0;
    mBeta  = 1.0;

    mA.reset( new lama::CSRSparseMatrix<double>() );

    lama::MatrixCreator::buildPoisson( *mA, mDimension, mStencilType, mSizeX, mSizeY, mSizeZ );

    dmemo::DistributionPtr rowDistribution = mA->getRowDistributionPtr();
    dmemo::DistributionPtr colDistribution = mA->getColDistributionPtr();

    mX.reset( new lama::DenseVector<double>() );
    mX->setSameValue( colDistribution, 1.0 );
    mY.reset( new lama::DenseVector<double>() );
    *mY = *mA * *mX;

    SCAI_LOG_INFO( logger, "created input set: A = " << *mA << ", X = " << *mX << ", Y = " << *mY );
}

benchmark::InputSet* PoissonInputSet::create( const std::string argument )
{
    // argment specfies kind of stencil, e.g. 2D5P_100_100

    return new PoissonInputSet( argument );
}

std::string PoissonInputSet::createValue()
{
    std::string id = "Poisson";
    return id;
}

const std::string& PoissonInputSet::getCreateId() const
{
    static std::string id = createValue();
    return id;
}

const std::string& PoissonInputSet::getArgument() const
{
    return mArgument;
}

}

}
