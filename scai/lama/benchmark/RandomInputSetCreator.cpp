/**
 * @file RandomInputSetCreator.cpp
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
 * @brief RandomInputSetCreator.cpp
 * @author Jiri Kraus
 * @date 06.05.2010
 * $Id$
 */

#include <scai/lama/benchmark/RandomInputSetCreator.hpp>

#include <scai/lama/matrix/CSRSparseMatrix.hpp>
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

using namespace scai;

// Random Creator

SCAI_LOG_DEF_LOGGER( RandomInputSetCreator::logger, "InputSetCreator.RandomInputSetCreator" );

RandomInputSetCreator::RandomInputSetCreator()
{
    SCAI_LOG_INFO( logger, "RandomInputSetCreator" );
}

RandomInputSetCreator::~RandomInputSetCreator()
{
    SCAI_LOG_INFO( logger, "~RandomInputSetCreator" );
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

LAMAInputSet* RandomInputSetCreator::create() const
{
    throw bf::BFError( "RandomInputSetCreator needs an argument of type size to create a matrix of size n" );
}

LAMAInputSet* RandomInputSetCreator::create( const std::string& arguments ) const
{
    SCAI_LOG_INFO( logger, "create " << arguments );

    std::istringstream input( arguments );

    IndexType size = 0;
    double fillingGrade = 1.0;

    input >> size >> fillingGrade;

    SCAI_LOG_INFO( logger, "create( " << arguments << ") : size = " << size << ", fillingGrade = " << fillingGrade );

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

    common::unique_ptr<lama::CSRSparseMatrix<double> > matrixA( new lama::CSRSparseMatrix<double>() );

    lama::MatrixCreator::buildRandom( *matrixA, size, fillingGrade );

    dmemo::DistributionPtr rowDistribution = matrixA->getRowDistributionPtr();
    dmemo::DistributionPtr colDistribution = matrixA->getColDistributionPtr();

    common::unique_ptr<lama::DenseVector<double> > vectorX( new lama::DenseVector<double>( colDistribution, 1.0 ) );
    common::unique_ptr<lama::DenseVector<double> > vectorY( new lama::DenseVector<double>( rowDistribution, 0.0 ) );

    SCAI_LOG_DEBUG( logger, "distributed vectors X, Y allocated: " << *vectorX << ", " << *vectorY );

    // Note: in previous versions Y was set when building Poisson matrices, now we do just mult

    *vectorY = ( *matrixA ) * ( *vectorX );

    SCAI_LOG_INFO( logger, "created input set: A = " << *matrixA << ", X = " << *vectorX << ", Y = " << *vectorY );

    return new LAMAInputSet( id(), 1.0, 1.0, vectorX, vectorY, matrixA );
}

LAMA_INPUTSET_REGISTRATION( RandomInputSetCreator );
