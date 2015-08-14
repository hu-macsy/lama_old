/**
 * @file SingleGridSetup.hpp
 *
 * @license
 * Copyright (c) 2009-2015
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
 * @brief SingleGridSetup.hpp
 * @author Jiri Kraus
 * @date 27.10.2011
 * @since 1.0.0
 */

// hpp
#include <scai/lama/solver/SingleGridSetup.hpp>

// others
#include <scai/lama/solver/SOR.hpp>
#include <scai/lama/solver/criteria/IterationCount.hpp>

// tracing
#include <scai/tracing.hpp>

namespace scai
{

namespace lama
{

SCAI_LOG_DEF_LOGGER( SingleGridSetup::logger, "AMGSetup.SingleGridSetup" )

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
    if( !mSolver )
    {
        SCAI_LOG_DEBUG( logger, "new sor" )
        SOR* sorSolver = new SOR( "10x SingleGridSetup SOR Solver" );

        CriterionPtr criterion( new IterationCount( 10 ) );

        sorSolver->setStoppingCriterion( criterion );

        mSolver.reset( sorSolver );
    }

    SCAI_LOG_DEBUG( logger, "mSolver->initialize" )
    mSolver->initialize( coefficients );

    SCAI_LOG_DEBUG( logger, "mIdentity.reset" )
    mIdentity.reset( coefficients.clone() );

    SCAI_LOG_DEBUG( logger, "before identity" )
    mIdentity->setIdentity( coefficients.getDistributionPtr() );
    SCAI_LOG_DEBUG( logger, "after identity" )

    SCAI_LOG_DEBUG( logger, "Identity matrix = " << *mIdentity )

    mSolutionVector.reset( Vector::createVector( coefficients.getValueType(), mIdentity->getDistributionPtr() ) );
    mRhsVector.reset( Vector::createVector( coefficients.getValueType(), mIdentity->getDistributionPtr() ) );
    mTmpResVector.reset( Vector::createVector( coefficients.getValueType(), mIdentity->getDistributionPtr() ) );
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

} /* end namespace lama */

} /* end namespace scai */
