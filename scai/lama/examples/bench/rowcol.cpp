/**
 * @file lama/examples/bench/rowcol.cpp
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
 * @brief Benchmark of matrix get/set row/col operations
 * @author Thomas Brandes
 * @date 12.10.2016
 */

#include <scai/lama.hpp>

#include <scai/tracing.hpp>
#include <scai/dmemo.hpp>
#include <scai/dmemo/BlockDistribution.hpp>
#include <scai/dmemo/NoDistribution.hpp>

#include <scai/common/Walltime.hpp>
#include <scai/common/Settings.hpp>

using namespace scai;
using namespace lama;

int main( int argc, const char* argv[] )
{
    scai::common::Settings::parseArgs( argc, const_cast<const char**>( argv ) );

    dmemo::CommunicatorPtr comm = dmemo::Communicator::getCommunicatorPtr();

    SCAI_REGION( "Main.BenchGetSet" )

    IndexType size = 1000;

    hmemo::ContextPtr ctx = hmemo::Context::getContextPtr();

    DenseMatrix<DefaultReal> mat( size, size );

    mat.setContextPtr( ctx );

    {
        SCAI_REGION( "Main.fillRandom" )
        MatrixCreator::fillRandom( mat, 0.1f );
    }

    dmemo::DistributionPtr dist( new dmemo::BlockDistribution( size, comm ) );
    dmemo::DistributionPtr rep( new dmemo::NoDistribution( size ) );

    mat.redistribute( dist, dist );

    DenseVector<DefaultReal> row( ctx );

    double tstart = common::Walltime::get();

    for ( IndexType i = 0; i < size; ++i )
    {
        mat.getRow( row, i );
        mat.setRow( row, i, common::BinaryOp::SUB );
    }

    double time = common::Walltime::get() - tstart;

    std::cout << "Time (set/get row) = " << time << " seconds" << std::endl;

    SCAI_ASSERT_LT( mat.maxNorm(), 0.001, "set/get row does not work correctly" );

    mat.allocate( size, size );

    {
        SCAI_REGION( "Main.fillRandom" )
        MatrixCreator::fillRandom( mat, 0.1f );
    }

    DenseVector<DefaultReal> col( ctx );

    mat.redistribute( dist, dist );

    tstart = common::Walltime::get();

    for ( IndexType j = 0; j < size; ++j )
    {
        mat.getColumn( col, j );
        mat.setColumn( col, j, common::BinaryOp::SUB );
    }

    time = common::Walltime::get() - tstart;

    std::cout << "Time (set/get column) = " << time << " seconds" << std::endl;

    SCAI_ASSERT_LT( mat.maxNorm(), 0.001, "set/get column does not work correctly" );
}

