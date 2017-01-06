/**
 * @file solver/SingleGridSetup.cpp
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
 * @brief Implementation of methods for the class SingleGridSetup.
 * @author Jiri Kraus
 * @date 27.10.2011
 */

// hpp
#include <scai/solver/SingleGridSetup.hpp>

// local library
#include <scai/solver/Jacobi.hpp>
#include <scai/solver/criteria/IterationCount.hpp>

// tracing
#include <scai/tracing.hpp>

namespace scai
{

namespace solver
{

SCAI_LOG_DEF_LOGGER( SingleGridSetup::logger, "AMGSetup.SingleGridSetup" )

using lama::Matrix;
using lama::Vector;
using lama::Scalar;

std::string SingleGridSetup::createValue()
{
    return "SingleGridSetup";
}

AMGSetup* SingleGridSetup::create()
{
    return new SingleGridSetup();
}

SingleGridSetup::SingleGridSetup()
{
    SCAI_LOG_DEBUG( logger, "SingleGridSetup" )
}

SingleGridSetup::~SingleGridSetup()
{
}

void SingleGridSetup::initialize( const Matrix& coefficients )
{
    SCAI_REGION( "initialize_SingleGridSetup" )
    SCAI_LOG_DEBUG( logger, "SingleGridSetup::initialize" )

    // set default solver
    if ( !mSolver )
    {
        SCAI_LOG_DEBUG( logger, "new Jacobi" )
        Jacobi* jacobiSolver = new Jacobi( "10x SingleGridSetup Jacobi Solver" );
        CriterionPtr criterion( new IterationCount( 10 ) );
        jacobiSolver->setStoppingCriterion( criterion );
        mSolver.reset( jacobiSolver );
    }

    SCAI_LOG_DEBUG( logger, "mSolver->initialize" )
    mSolver->initialize( coefficients );
    SCAI_LOG_DEBUG( logger, "mIdentity.reset" )
    mIdentity.reset( coefficients.newMatrix() );
    SCAI_LOG_DEBUG( logger, "before identity" )
    mIdentity->setIdentity( coefficients.getRowDistributionPtr() );
    SCAI_LOG_DEBUG( logger, "after identity" )
    SCAI_LOG_DEBUG( logger, "Identity matrix = " << *mIdentity )
    mSolutionVector.reset( coefficients.newDenseVector() );
    mRhsVector.reset( coefficients.newDenseVector() );
    mTmpResVector.reset( coefficients.newDenseVector() );
}

Solver& SingleGridSetup::getCoarseLevelSolver()
{
    return *mSolver;
}

unsigned int SingleGridSetup::getNumLevels()
{
    return 2;
}

Solver& SingleGridSetup::getSmoother( const unsigned int )
{
    return *mSolver;
}

const Matrix& SingleGridSetup::getGalerkin( const unsigned int )
{
    return mSolver->getCoefficients();
}

const Matrix& SingleGridSetup::getRestriction( const unsigned int )
{
    return *mIdentity;
}

const Matrix& SingleGridSetup::getInterpolation( const unsigned int )
{
    return *mIdentity;
}

Vector& SingleGridSetup::getSolutionVector( const unsigned int )
{
    return *mSolutionVector;
}

Vector& SingleGridSetup::getRhsVector( const unsigned int )
{
    return *mRhsVector;
}

Vector& SingleGridSetup::getTmpResVector( const unsigned int )
{
    return *mTmpResVector;
}

std::string SingleGridSetup::getCouplingPredicateInfo() const
{
    return "No coupling predicate.";
}

std::string SingleGridSetup::getColoringInfo() const
{
    return "No coloring.";
}

std::string SingleGridSetup::getInterpolationInfo() const
{
    return "Identity.";
}

std::string SingleGridSetup::getSmootherInfo() const
{
    return mSolver->getId();
}

std::string SingleGridSetup::getCoarseLevelSolverInfo() const
{
    return mSolver->getId();
}

void SingleGridSetup::setCoarseLevelSolver( SolverPtr solver )
{
    mSolver = solver;
}

void SingleGridSetup::setSmoother( SolverPtr solver )
{
    mSolver = solver;
}

void SingleGridSetup::writeAt( std::ostream& stream ) const
{
    stream << "SingleGridSetup";
}

} /* end namespace solver */

} /* end namespace scai */
