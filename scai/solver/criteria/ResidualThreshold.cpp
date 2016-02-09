/**
 * @file ResidualThreshold.cpp
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
 * @brief ResidualThreshold.cpp
 * @author Kai Buschulte
 * @date 21.07.2011
 * @since 1.0.0
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
    SCAI_ASSERT( normResult >= 0.0, "A norm should be always positive but is " << normResult );

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
