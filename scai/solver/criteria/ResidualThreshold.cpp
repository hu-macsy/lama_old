/**
 * @file ResidualThreshold.cpp
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
 * @endlicense
 *
 * @brief ResidualThreshold.cpp
 * @author Kai Buschulte
 * @date 21.07.2011
 */

// hpp
#include <scai/solver/criteria/ResidualThreshold.hpp>

// local library
#include <scai/solver/IterativeSolver.hpp>

// internal libraries
#include <scai/lama/norm/L2Norm.hpp> //default

namespace scai
{

namespace solver
{

ResidualThreshold::ResidualThreshold()
    : Criterion(), mNorm( lama::NormPtr( new lama::L2Norm() ) ), mCheckMode( Relative ), mPrecision( lama::Scalar( 1e-5 ) ), mFirstNormResult(
          -1.0 )
{
}

ResidualThreshold::ResidualThreshold( const lama::NormPtr norm )
    : Criterion(), mNorm( norm ), mCheckMode( Relative ), mPrecision( lama::Scalar( 1e-5 ) ), mFirstNormResult(
          -1.0 )
{
}

ResidualThreshold::ResidualThreshold( const lama::NormPtr norm, lama::Scalar precision, ResidualThresholdCheckMode checkMode )
    : Criterion(), mNorm( norm ), mCheckMode( checkMode ), mPrecision( precision ), mFirstNormResult( -1.0 )

{
}

ResidualThreshold::ResidualThreshold( const ResidualThreshold& other )
    : Criterion(), mNorm( other.mNorm ), mCheckMode( other.mCheckMode ), mPrecision( other.mPrecision ), mFirstNormResult(
          other.mFirstNormResult )
{
}

ResidualThreshold::~ResidualThreshold()
{
}

Criterion* ResidualThreshold::copy() const
{
    return new ResidualThreshold( *this );
}

inline bool ResidualThreshold::isSatisfied( const IterativeSolver& solver )
{
    lama::Scalar normResult = ( *mNorm )( solver.getResidual() );
    SCAI_ASSERT( normResult > 0.0 || normResult == 0.0, "A norm should be always positive but is " << normResult );

    switch( mCheckMode )
    {
        case ResidualThreshold::Absolute:
            SCAI_LOG_DEBUG( logger,
                            "Absolute residual in iteration " << solver.getIterationCount() << " is " << normResult << " should become smaller than precision " << mPrecision )
            ;
            return normResult < mPrecision;

        case ResidualThreshold::Relative:
        {
            if( mFirstNormResult == -1.0 )
            {
                //TODO define member variable for solver with getInitialResidual function
                mFirstNormResult = normResult;
            }

            SCAI_LOG_DEBUG( logger,
                            "Relative residual in iteration " << solver.getIterationCount() << " is " << normResult << " divided by firstNormResult " << mFirstNormResult << " is " << normResult/mFirstNormResult << " should become smaller than precision " << mPrecision );
            return ( normResult / mFirstNormResult ) < mPrecision;
        }

        case ResidualThreshold::Divergence:
        {
            if( mFirstNormResult == -1.0 )
            {
                //TODO define member variable for solver with getInitialResidual function
                mFirstNormResult = normResult;
            }

            SCAI_LOG_DEBUG( logger,
                            "Relative residual in iteration " << solver.getIterationCount() << " is " << normResult << " divided by firstNormResult " << mFirstNormResult << " is " << normResult/mFirstNormResult << " should become larger than precision " << mPrecision );
            return ( normResult / mFirstNormResult ) > mPrecision;
        }
    }

    return false;
}

ResidualThreshold::ResidualThresholdCheckMode ResidualThreshold::getCheckMode() const
{
    return mCheckMode;
}

lama::Scalar ResidualThreshold::getFirstNormResult() const
{
    return mFirstNormResult;
}

const lama::NormPtr ResidualThreshold::getNorm() const
{
    return mNorm;
}

lama::Scalar ResidualThreshold::getPrecision() const
{
    return mPrecision;
}

void ResidualThreshold::setCheckMode( ResidualThresholdCheckMode checkMode )
{
    mCheckMode = checkMode;
}

void ResidualThreshold::setFirstNormResult( lama::Scalar firstNormResult )
{
    mFirstNormResult = firstNormResult;
}

void ResidualThreshold::setPrecision( lama::Scalar precision )
{
    mPrecision = precision;
}

void ResidualThreshold::writeAt( std::ostream& stream ) const
{
    stream << "ResThr<" << getPrecision() << ", " << mCheckMode << ">";
}

} /* end namespace solver */

} /* end namespace scai */
