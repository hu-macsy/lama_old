/**
 * @file RandomInputSet.cpp
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
 * @brief RandomInputSet.cpp
 * @author Thomas Brandes
 * @date 15.09.2017
 */

#include <scai/benchmark/Parser.hpp>

#include <scai/lama/benchmark/RandomInputSet.hpp>
#include <scai/lama/expression/MatrixVectorExpressions.hpp>
#include <scai/lama/matutils/MatrixCreator.hpp>

namespace scai
{

namespace lama
{

RandomInputSet::RandomInputSet( const std::string arguments ) : LAMAInputSet()
{
    SCAI_LOG_INFO( logger, "create " << arguments );

    std::vector<std::string> args;

    tokenize( args, arguments, " ,:" );
 
    SCAI_ASSERT_EQ_ERROR( args.size(), 2, "arguments: " << arguments );

    IndexType size = std::stoi( args[0] );
    double fillingGrade = std::stod( args[1] );

    SCAI_LOG_INFO( logger, "create( " << arguments << ") : size = " << size << ", fillingGrade = " << fillingGrade );

    if ( fillingGrade <= 0.0 || fillingGrade > 1.0 ) 
    {
        throw benchmark::BFException(
            "RandomInputSetCreator creator needs valid filling grade > 0.0 and <= 1.0 to create a InputSet." );
    }

    mA.reset( new lama::CSRSparseMatrix<double>() );

    SCAI_LOG_INFO( logger, "mA = " << *mA )

    lama::MatrixCreator::buildRandom( *mA, size, fillingGrade );

    SCAI_LOG_INFO( logger, "random: mA = " << *mA )

    dmemo::DistributionPtr colDistribution = mA->getColDistributionPtr();

    SCAI_LOG_INFO( logger, "col dist = " << *colDistribution )

    mX.reset( new lama::DenseVector<double>( colDistribution, 1.0 ) );

    SCAI_LOG_INFO( logger, "mX = " << *mX )

    mY.reset( new lama::DenseVector<double>( *mA * *mX ) );

    SCAI_LOG_INFO( logger, "created input set: A = " << *mA << ", X = " << *mX << ", Y = " << *mY );
}

benchmark::InputSet* RandomInputSet::create( const std::string argument )
{
    // argment is take as correspoding file name

    return new RandomInputSet( argument );
}

std::string RandomInputSet::createValue()
{
    std::string id = "Random";
    return id;
}

}

}
