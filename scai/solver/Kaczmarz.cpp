/**
 * @file Kaczmarz.cpp
 *
 * @license
 * Copyright (c) 2009-2016
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
 * @brief Implementation of methods for the Kaczmarz solver.
 * @author Thomas Brandes
 * @date 24.08.2011
 */

// hpp
#include <scai/solver/Kaczmarz.hpp>

// internal scai libraries
#include <scai/lama/DenseVector.hpp>

#include <scai/lama/expression/VectorExpressions.hpp>
#include <scai/lama/expression/MatrixVectorExpressions.hpp>

#include <scai/tracing.hpp>

#include <scai/common/SCAITypes.hpp>

namespace scai
{

namespace solver
{

SCAI_LOG_DEF_LOGGER( Kaczmarz::logger, "Solver.IterativeSolver.Kaczmarz" )

using lama::Matrix;
using lama::Vector;
using lama::Scalar;

Kaczmarz::Kaczmarz( const std::string& id )
    : IterativeSolver( id )
{
}

Kaczmarz::Kaczmarz( const std::string& id, LoggerPtr logger )
    : IterativeSolver( id, logger )
{
}

Kaczmarz::Kaczmarz( const Kaczmarz& other )
    : IterativeSolver( other )
{
}

Kaczmarz::~Kaczmarz()
{
}

Kaczmarz::KaczmarzRuntime::KaczmarzRuntime()

    : IterativeSolverRuntime()
{
}

Kaczmarz::KaczmarzRuntime::~KaczmarzRuntime()
{
}

void Kaczmarz::initialize( const Matrix& coefficients )
{
    SCAI_REGION( "Solver.Kaczmarz.initialize" )
    IterativeSolver::initialize( coefficients );
    KaczmarzRuntime& runtime = getRuntime();
    runtime.mRow.reset( coefficients.newDenseVector() );
}

void Kaczmarz::iterate()
{
    SCAI_REGION( "Solver.Kaczmarz.iterate" )

    KaczmarzRuntime& runtime = getRuntime();

    IndexType iter = this->getIterationCount();

    if ( iter == 0 )
    {
        this->getResidual();
    }

    const Matrix& A = *runtime.mCoefficients;
    const Vector& b = *runtime.mRhs;
    Vector& x = *runtime.mSolution;
    Vector& z = *runtime.mRow;

    SCAI_LOG_INFO( logger, "Iteration " << iter )

    // Kaczmarz implementation start

    for ( IndexType iRow = 0; iRow < A.getNumRows(); ++iRow )
    {
        SCAI_LOG_DEBUG( logger, "Update irow " << iRow )

        A.getRow( z, iRow );

        z.redistribute( x.getDistributionPtr() );

        Scalar bi = b( iRow );

        Scalar p = ( bi - z.dotProduct( x ) ) / z.dotProduct( z );

        x = x + p * z;
    }


    // Kaczmarz implementation end

    mKaczmarzRuntime.mSolution.setDirty( true );
}

SolverPtr Kaczmarz::copy()
{
    return SolverPtr( new Kaczmarz( *this ) );
}

Kaczmarz::KaczmarzRuntime& Kaczmarz::getRuntime()
{
    return mKaczmarzRuntime;
}

const Kaczmarz::KaczmarzRuntime& Kaczmarz::getConstRuntime() const
{
    return mKaczmarzRuntime;
}

std::string Kaczmarz::createValue()
{
    return "Kaczmarz";
}

Solver* Kaczmarz::create( const std::string name )
{
    return new Kaczmarz( name );
}

void Kaczmarz::writeAt( std::ostream& stream ) const
{
    stream << "Kaczmarz ( id = " << mId << ", #iter = " << getConstRuntime().mIterations << " )";
}

} /* end namespace solver */

} /* end namespace scai */
