/**
 * @file kaczmarz.cpp
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
 * @brief Kaczmarz method with sparse vector
 * @author Thomas Brandes
 * @date 23.02.2017
 */

#include <scai/lama.hpp>

#include <scai/solver/logger/CommonLogger.hpp>
#include <scai/solver/criteria/IterationCount.hpp>
#include <scai/solver/criteria/ResidualThreshold.hpp>
#include <scai/solver/CG.hpp>

// _Matrix & vector related includes
#include <scai/lama/DenseVector.hpp>
#include <scai/lama/SparseVector.hpp>
#include <scai/lama/matrix/CSRSparseMatrix.hpp>
#include <scai/lama/matrix/DenseMatrix.hpp>
#include <scai/lama/storage/CSRStorage.hpp>
#include <scai/dmemo/BlockDistribution.hpp>

// import common 
#include <scai/common/Walltime.hpp>
#include <scai/common/Settings.hpp>

#include <iostream>
#include <stdlib.h>

using namespace scai;
using namespace lama;
using namespace solver;

typedef DefaultReal ValueType;

int main( int argc, const char* argv[] )
{
    // relevant SCAI arguments: 
    //   SCAI_CONTEXT = ...    set default context
    //   SCAI_DEVICE  = ...    set default device

    common::Settings::parseArgs( argc, argv );

    std::string filename = argv[1];

    auto A = read<CSRSparseMatrix<ValueType>>( filename );
  
    DenseVector<ValueType> diag;

    ValueType theta = 0.0;

    A.getDiagonal( diag );
    diag -= theta;
    A.setDiagonal( diag );

    auto x = fill<DenseVector<ValueType>>( A.getRowDistributionPtr(), 1 );

    DenseVector<ValueType> q;   // tmp vector, allocated data will be reused in loop
    DenseVector<ValueType> y;   // tmp vector, allocated data will be reused in loop

    LoggerPtr slogger( new CommonLogger( "CGLogger:", LogLevel::convergenceHistory, LoggerWriteBehaviour::toConsoleOnly ) );
    CG<ValueType> solver( "CG", slogger );
    solver.initialize( A );

    NormPtr<ValueType> norm( new L1Norm<ValueType>() );
    ValueType eps = 1e-6;

    CriterionPtr<ValueType> rt( new ResidualThreshold<ValueType>( norm, eps, ResidualCheck::Absolute ) );
    CriterionPtr<ValueType> it( new IterationCount<ValueType>( 100 ) );

    // stop if iteration count reached OR residual threshold is reached

    CriterionPtr<ValueType> criterion( new Criterion<ValueType>( it, rt, BooleanOp::OR ) );
    solver.setStoppingCriterion( criterion );

    for ( int k = 0; k < 10; ++k )
    {
        ValueType norm = x.l2Norm();
        q = x / norm;
        solver.solve( x, q );
        y = A * x;
        ValueType lambda = x.dotProduct ( y ) / x.dotProduct( x );
        std::cout << "lambda = " << lambda << std::endl;
    }

    x.writeToFile( "eigenvector.mtx" );
}
