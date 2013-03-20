/**
 * @file LAMAInputSetCreators.cpp
 *
 * @license
 * Copyright (c) 2011
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
 * @brief LAMAInputSetCreators.cpp
 * @author Jiri Kraus
 * @date 06.05.2010
 * $Id$
 */

#include <bench/LAMAInputSetCreators.hpp>

#include <framework/src/benchmark_framework.hpp>

#include <lama/matrix/CSRSparseMatrix.hpp>

#include <lama/HostWriteAccess.hpp>
#include <lama/CommunicatorFactory.hpp>
#include <lama/distribution/BlockDistribution.hpp>
#include <lama/distribution/GeneralDistribution.hpp>
#include <lama/matutils/MatrixCreator.hpp>
#include <lama/expression/MatrixVectorExpressions.hpp>

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
    } while( found != std::string::npos );

    return args;
}
}

LAMA_LOG_DEF_LOGGER( PoissonInputSetCreator::logger, "InputSetCreator.PoissonInputSetCreator" );

PoissonInputSetCreator::PoissonInputSetCreator()
{
    LAMA_LOG_INFO( logger, "PoissonInputSetCreator" );
}

PoissonInputSetCreator::~PoissonInputSetCreator()
{
    LAMA_LOG_INFO( logger, "~PoissonInputSetCreator" );
}

std::auto_ptr<PoissonInputSetCreator::InputSetType> PoissonInputSetCreator::create() const
{
    std::ostringstream message;
    message << "PoissonInputSetCreator::create: This InputSetCreator requires the following arguments:" << std::endl
            << "    aDbP_dimX_dimY_dimZ," << std::endl
            << "where a is in (1, 2, 3) and b is in (3, 5, 7, 9, 19, 27) and dimY "
            << "and dimZ are needed in case of unequality to one, only." << std::endl;
    throw bf::BFException( message.str() );
}

std::auto_ptr<PoissonInputSetCreator::InputSetType> PoissonInputSetCreator::createSet( const std::string& params )
{
    LAMA_LOG_INFO( logger, "create, params = " << params );

    // 2D5P_10000_10
    std::vector<std::string> args = split( params, '_' );

    if( args.size() < 2 || args.size() > 4 )
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

    if( args[0].find( 'D', 0 ) == std::string::npos || args[0].find( 'P', 0 ) == std::string::npos )
    {
        std::ostringstream message;
        message << "PoissonInputSetCreator::create( ): Expected 'D' and 'P' within " << args[0]
                << " but at least one was not found!" << std::endl;
        throw bf::BFException( message.str() );
    }

    vals = split( args[0], 'D' );

    if( vals.size() != 2 )
    {
        std::ostringstream message;
        message << "PoissonInputSetCreator::create( ): Illegal format for " << args[0] << std::endl
                << "Should have been aDbP, with a=1,2,3 and b=3,5,7,9,19 or 27." << std::endl;
        throw bf::BFException( message.str() );
    }

    vals[1] = split( vals[1], 'P' )[0];
    std::istringstream convert( vals[0] );
    convert >> dimension;

    LAMA_LOG_DEBUG( logger, "Stencil dimension = " << dimension );

    if( convert.fail() )
    {
        std::ostringstream message;
        message << "PoissonInputSetCreator::create( ): Unable to convert '" << convert.str()
                << "' (should have been dimension)." << std::endl;
        throw bf::BFException( message.str() );
    }

    convert.clear();
    convert.str( vals[1] );
    convert >> stencilType;

    if( convert.fail() )
    {
        std::ostringstream message;
        message << "PoissonInputSetCreator::create( ): Unable to convert '" << convert.str()
                << "' (should have been stencilType)." << std::endl;
        throw bf::BFException( message.str() );
    }

    convert.clear();

    if( dimension != (IndexType) args.size() - 1 )
    {
        std::ostringstream message;
        message << "PoissonInputSetCreator::create( ): '" << params << "' said " << "to have " << dimension
                << " dimensions, but had " << ( args.size() - 1 ) << " dimension parameters." << std::endl;
        throw bf::BFException( message.str() );
    }

    switch( dimension )
    {
    case 3:
        convert.str( args[3] );
        convert >> std::dec >> dimZ;

        if( convert.fail() )
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

        if( convert.fail() )
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

        if( convert.fail() )
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

    LAMA_LOG_INFO( logger, "dimX = " << dimX << ", dimY = " << dimY << ", dimZ = " << dimZ );

    bool exception = false;

    if( dimZ == 1 && dimY == 1 )
    {
        exception = ( stencilType != 3 );
    }
    else if( dimZ == 1 )
    {
        exception = ( stencilType != 5 && stencilType != 9 );
    }
    else if( dimZ >= 1 && dimY >= 1 )
    {
        exception = ( stencilType != 7 && stencilType != 19 && stencilType != 27 );
    }
    else
    {
        exception = true;
    }

    if( exception )
    {
        std::ostringstream message;
        message << "PoissonInputSetCreator::create( ): Wrong values for parameters:" << std::endl << "Valid values are:"
                << std::endl << " dimY  | dimZ  | stencilType " << std::endl << "-------+-------+-------------"
                << std::endl << "  == 1 |  == 1 |           3 " << std::endl << "  >  1 |  == 1 |        5, 9 "
                << std::endl << "  >= 1 |  >= 1 |   7, 19, 27 " << std::endl;
        throw bf::BFException( message.str() );
    }

    std::auto_ptr<lama::CSRSparseMatrix<double> > matrixA( new lama::CSRSparseMatrix<double>() );

    lama::MatrixCreator<double>::buildPoisson( *matrixA, dimension, stencilType, dimX, dimY, dimZ );

    lama::DistributionPtr rowDistribution = matrixA->getDistributionPtr();
    lama::DistributionPtr colDistribution = matrixA->getColDistributionPtr();

    std::auto_ptr<lama::DenseVector<double> > vectorX( new lama::DenseVector<double>( colDistribution, 1.0 ) );
    std::auto_ptr<lama::DenseVector<double> > vectorY( new lama::DenseVector<double>( rowDistribution, 0.0 ) );

    LAMA_LOG_DEBUG( logger, "distributed vectors X, Y allocated: " << *vectorX << ", " << *vectorY );

    // Note: in previous versions Y was set when building Poisson matrices, now we do just mult

    *vectorY = ( *matrixA ) * ( *vectorX );

    LAMA_LOG_INFO( logger, "created input set: A = " << *matrixA << ", X = " << *vectorX << ", Y = " << *vectorY );

    std::auto_ptr<InputSetType>

    iSet( new InputSetType( id(), 1.0, 1.0, vectorX, vectorY, matrixA ) );

    return iSet;
}

std::auto_ptr<PoissonInputSetCreator::InputSetType> PoissonInputSetCreator::create( const std::string& params ) const
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

// EmptyDiagonalCreator

LAMA_LOG_DEF_LOGGER( EmptyDiagonalInputSetCreator::logger, "InputSetCreator.EmptyDiagonalInputSetCreator" );

EmptyDiagonalInputSetCreator::EmptyDiagonalInputSetCreator()
{
    LAMA_LOG_INFO( logger, "EmptyDiagonalInputSetCreator" );
}

EmptyDiagonalInputSetCreator::~EmptyDiagonalInputSetCreator()
{
    LAMA_LOG_INFO( logger, "~EmptyDiagonalInputSetCreator" );
}

const std::string& EmptyDiagonalInputSetCreator::id()
{
    static const std::string id = "EmptyDiagonal";
    return id;
}

const std::string& EmptyDiagonalInputSetCreator::getId() const
{
    return id();
}

std::auto_ptr<LAMAInputSet> EmptyDiagonalInputSetCreator::create() const
{
    throw bf::BFError( "EmptyDiagonalInputSetCreator needs an argument of type size to create a matrix of size n" );
}

namespace
{
inline lama::IndexType sum( lama::IndexType n )
{
    return ( n * ( n + 1 ) ) >> 1;
}
}

std::auto_ptr<LAMAInputSet> EmptyDiagonalInputSetCreator::create( const std::string& arguments ) const
{
    LAMA_LOG_INFO( logger, "create " << arguments );

    std::istringstream convert( arguments );
    IndexType size;
    convert >> std::dec >> size;

    if( convert.fail() )
    {
        throw new bf::BFError(
            "Invalid argument '" + arguments + "'. Need just one argument for the size of the matrix (e.g. "
            + id() + "( 2048 )." );
    }

    if( size <= 0 )
    {
        throw new bf::BFError(
            "Argument " + arguments + " was an invalid integer. Size needs to be greater than zero!" );
    }

    ValueType* const values = new ValueType[size * size];
    ValueType* const rhs = new ValueType[size];

    for( IndexType i = 0; i < size; ++i )
    {
        for( IndexType j = 0; j < size; ++j )
        {
            values[i * size + j] = std::fabs( static_cast<ValueType>( i - j ) );
        }

        rhs[i] = sum( size - i - 1 ) + sum( i );
    }

    std::auto_ptr<lama::DenseVector<double> > vecX( new lama::DenseVector<double>( size, 1.0 ) );
    std::auto_ptr<lama::DenseVector<double> > vecY( new lama::DenseVector<double>( size, rhs ) );

    std::auto_ptr<lama::CSRSparseMatrix<double> > matrix( new lama::CSRSparseMatrix<double>() );
    matrix->setRawDenseData( size, size, values );
    delete[] values;
    delete[] rhs;
    std::auto_ptr<LAMAInputSet> iSet( new LAMAInputSet( getId(), 1.0, 1.0, vecX, vecY, matrix ) );
    return iSet;
}

LAMA_INPUTSET_REGISTRATION( EmptyDiagonalInputSetCreator );

// Random Creator

LAMA_LOG_DEF_LOGGER( RandomInputSetCreator::logger, "InputSetCreator.RandomInputSetCreator" );

RandomInputSetCreator::RandomInputSetCreator()
{
    LAMA_LOG_INFO( logger, "RandomInputSetCreator" );
}

RandomInputSetCreator::~RandomInputSetCreator()
{
    LAMA_LOG_INFO( logger, "~RandomInputSetCreator" );
}

const std::string& RandomInputSetCreator::id()
{
    static const std::string id = "Random";
    return id;
}

const std::string& RandomInputSetCreator::getId() const
{
    return id();
}

std::auto_ptr<LAMAInputSet> RandomInputSetCreator::create() const
{
    throw bf::BFError( "RandomInputSetCreator needs an argument of type size to create a matrix of size n" );
}

std::auto_ptr<LAMAInputSet> RandomInputSetCreator::create( const std::string& arguments ) const
{
    LAMA_LOG_INFO( logger, "create " << arguments );

    std::istringstream input( arguments );

    IndexType size = 0;
    double fillingGrade = 1.0;

    input >> size >> fillingGrade;

    LAMA_LOG_INFO( logger, "create( " << arguments << ") : size = " << size << ", fillingGrade = " << fillingGrade );

    if( size <= 0 )
    {
        throw bf::BFException( "RandomInputSetCreator creator needs a size to create a InputSet." );
    };

    if( input.good() && ( fillingGrade <= 0.0 || fillingGrade > 1.0 ) )
    {
        throw bf::BFException(
            "RandomInputSetCreator creator needs valid filling grade > 0.0 and <= 1.0 to create a InputSet." );
    }

    if( input.fail() )
    {
        fillingGrade = 1.0;
    }

    std::auto_ptr<lama::CSRSparseMatrix<double> > matrixA( new lama::CSRSparseMatrix<double>() );

    lama::MatrixCreator<double>::buildRandom( *matrixA, size, fillingGrade );

    lama::DistributionPtr rowDistribution = matrixA->getDistributionPtr();
    lama::DistributionPtr colDistribution = matrixA->getColDistributionPtr();

    std::auto_ptr<lama::DenseVector<double> > vectorX( new lama::DenseVector<double>( colDistribution, 1.0 ) );
    std::auto_ptr<lama::DenseVector<double> > vectorY( new lama::DenseVector<double>( rowDistribution, 0.0 ) );

    LAMA_LOG_DEBUG( logger, "distributed vectors X, Y allocated: " << *vectorX << ", " << *vectorY );

    // Note: in previous versions Y was set when building Poisson matrices, now we do just mult

    *vectorY = ( *matrixA ) * ( *vectorX );

    LAMA_LOG_INFO( logger, "created input set: A = " << *matrixA << ", X = " << *vectorX << ", Y = " << *vectorY );

    std::auto_ptr<LAMAInputSet> iSet( new LAMAInputSet( id(), 1.0, 1.0, vectorX, vectorY, matrixA ) );

    return iSet;
}

LAMA_INPUTSET_REGISTRATION( RandomInputSetCreator );
