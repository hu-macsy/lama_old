/**
 * @file RandomInputSetCreator.cpp
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
 * @brief RandomInputSetCreator.cpp
 * @author Jiri Kraus
 * @date 06.05.2010
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

namespace scai
{

namespace lama
{

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
    throw benchmark::BFError( "RandomInputSetCreator needs an argument of type size to create a matrix of size n" );
}

LAMAInputSet* RandomInputSetCreator::create( const std::string& arguments ) const
{
    SCAI_LOG_INFO( logger, "create " << arguments );

    std::istringstream input( arguments );

    IndexType size = 0;
    double fillingGrade = 1.0;

    input >> size >> fillingGrade;

    SCAI_LOG_INFO( logger, "create( " << arguments << ") : size = " << size << ", fillingGrade = " << fillingGrade );

    if ( size <= 0 )
    {
        throw benchmark::BFException( "RandomInputSetCreator creator needs a size to create a InputSet." );
    };

    if ( input.good() && ( fillingGrade <= 0.0 || fillingGrade > 1.0 ) )
    {
        throw benchmark::BFException(
            "RandomInputSetCreator creator needs valid filling grade > 0.0 and <= 1.0 to create a InputSet." );
    }

    if ( input.fail() )
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

}

}

