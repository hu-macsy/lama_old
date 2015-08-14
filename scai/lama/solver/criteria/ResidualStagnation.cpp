/**
 * @file ResidualStagnation.hpp
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
 * @brief ResidualStagnation.hpp
 * @author Kai Buschulte
 * @date 25.07.2011
 * @since 1.0.0
 */

// hpp
#include <scai/lama/solver/criteria/ResidualStagnation.hpp>

// others
#include <scai/lama/solver/IterativeSolver.hpp>

#include <iostream>
#include <algorithm>

namespace lama
{

ResidualStagnation::ResidualStagnation( NormPtr norm )
    : Criterion(), mNorm( norm ), mLookback( 1 ), //TODO check default value
      mLastResidualNorms( 1 ), mNextEntry( 0 ), mEntriesReady( false ), mPrecision( 0.1 ) //TODO check default value
{
}

ResidualStagnation::ResidualStagnation( NormPtr norm, IndexType lookback, Scalar precision )
    : Criterion(), mNorm( norm ), mLookback( lookback ), mLastResidualNorms( lookback ), mNextEntry( 0 ), mEntriesReady(
          false ), mPrecision( precision )
{
}

ResidualStagnation::ResidualStagnation( const ResidualStagnation &other )
    : Criterion(), mNorm( other.mNorm ), mLookback( other.mLookback ), mLastResidualNorms(
          other.mLastResidualNorms ), mNextEntry( other.mNextEntry ), mEntriesReady(
          other.mEntriesReady ), mPrecision( other.mPrecision )
{
}

ResidualStagnation::~ResidualStagnation()
{
}

Criterion* ResidualStagnation::copy() const
{
    return new ResidualStagnation( *this );
}

bool ResidualStagnation::isSatisfied( const IterativeSolver& solver )
{
    mLastResidualNorms[mNextEntry] = ( *mNorm )( solver.getResidual() ).getValue<long double>();

    mNextEntry = ( mNextEntry + 1 ) % mLookback;

    if( mNextEntry == 0 )
    {
        mEntriesReady = true;
    }

    if( mEntriesReady )
    {
        Scalar min = *std::min_element( mLastResidualNorms.begin(), mLastResidualNorms.end() );
        Scalar max = *std::max_element( mLastResidualNorms.begin(), mLastResidualNorms.end() );

        min = std::max( std::numeric_limits<Scalar>::min(), min );
        //std::cout<< " max ="<<max<<"       min = "<<min<<"      max/min = "<<max/min<<"1+p = "<<(1.0+mPrecision)<<std::endl;
        mEntriesReady = false;
        return ( ( max / min ) < ( 1.0 + mPrecision ) );
    }

    return false;
}

int ResidualStagnation::getLookback() const
{
    return mLookback;
}

const NormPtr ResidualStagnation::getNorm() const
{
    return mNorm;
}

const Scalar ResidualStagnation::getPrecision() const
{
    return mPrecision;
}

void ResidualStagnation::setLookback( IndexType lookback )
{
    mLookback = lookback;
}

void ResidualStagnation::setPrecision( const Scalar precision )
{
    mPrecision = precision;
}

void ResidualStagnation::writeAt( std::ostream& stream ) const
{
    stream << "ResStgn<" << getPrecision() << ", " << getLookback() << ">";
}

} // namespace lama
