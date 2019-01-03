/**
 * @file gather.cpp
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
 * @brief Example of global gathering of data.
 * @author Thomas Brandes
 * @date 02.01.2019
 */

#include <scai/lama.hpp>

#include <scai/lama/DenseVector.hpp>

#include <scai/dmemo/BlockDistribution.hpp>
#include <scai/dmemo/CyclicDistribution.hpp>
#include <scai/dmemo/GlobalAddressingPlan.hpp>
#include <scai/lama/DenseVector.hpp>
#include <scai/utilskernel/HArrayUtils.hpp>
#include <scai/common/Walltime.hpp>

#include <iostream>
#include <stdlib.h>

using namespace scai;
using namespace lama;
using namespace hmemo;
using common::Walltime;

typedef DefaultReal ValueType;

int main()
{
    const IndexType N = 100000;
    const IndexType M = 30000;

    auto sourceDist = dmemo::blockDistribution( N );
    auto targetDist = dmemo::cyclicDistribution( M, 2 );

    auto fillSource = []( IndexType k ) { return ValueType( k ) / ValueType( N ); };

    auto indexes = linearDenseVector<IndexType>( targetDist, 2, 1 );
    auto source  = fillDenseVector<ValueType>( sourceDist, fillSource );
    
    DenseVector<ValueType> target( targetDist, ValueType( 0 ) );
    target.gather( source, indexes, common::BinaryOp::COPY );

    SCAI_ASSERT_EQ_ERROR( target.getDistribution(), indexes.getDistribution(), "serious error" );

    const IndexType localM = targetDist->getLocalSize();

    HArray<IndexType> csrIA;
    utilskernel::HArrayUtils::setSequence( csrIA, 0, 1, localM + 1 );
    HArray<IndexType> csrJA( indexes.getLocalValues() );
    HArray<ValueType> csrValues( localM, ValueType( 1 ) );

    CSRStorage<ValueType> csrLocal( localM, N, std::move( csrIA ), std::move( csrJA ), std::move( csrValues ) );

    CSRSparseMatrix<ValueType> gatherMatrix( targetDist, std::move( csrLocal ) );
    gatherMatrix.redistribute( targetDist, sourceDist );

    DenseVector<ValueType> target1;
    target1 = gatherMatrix * source;

    auto maxDiff = target.maxDiffNorm( target1 );

    std::cout << "maxDiff = " << maxDiff << std::endl;

    auto plan = dmemo::globalAddressingPlan( *sourceDist, indexes.getLocalValues() );
    auto target2 = fillDenseVector<ValueType>( targetDist, 0 );
    plan.gather( target2.getLocalValues(), source.getLocalValues(), common::BinaryOp::COPY );

    auto maxDiff2 = target.maxDiffNorm( target2 );

    std::cout << "maxDiff = " << maxDiff2 << std::endl;

    const IndexType NITER = 1000;

    double timeV = Walltime::get();

    for ( IndexType k = 0; k < NITER; ++k )
    {
        target.gather( source, indexes, common::BinaryOp::COPY );
    }

    timeV = Walltime::get() - timeV;

    std::cout << "Time gatherVector = " << timeV << std::endl;

    double timeP = Walltime::get();

    for ( IndexType k = 0; k < NITER; ++k )
    {
        plan.gather( target.getLocalValues(), source.getLocalValues(), common::BinaryOp::COPY );
    }

    timeP = Walltime::get() - timeP;

    std::cout << "Time gatherVector (plan) = " << timeP << std::endl;

    double timeM = Walltime::get();

    for ( IndexType k = 0; k < NITER; ++k )
    {
         target1 = gatherMatrix * source;
    }

    timeM = Walltime::get() - timeM;
    std::cout << "Time gatherMatrix = " << timeM << std::endl;
}
