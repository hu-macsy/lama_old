/**
 * @file solver/criteria/ResidualStagnation.cpp
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
 * @brief ResidualStagnation.hpp
 * @author Kai Buschulte
 * @date 25.07.2011
 */

// hpp
#include <scai/solver/criteria/ResidualStagnation.hpp>

// local library
#include <scai/solver/IterativeSolver.hpp>

// std
#include <iostream>
#include <algorithm>

namespace scai
{

namespace solver
{

ResidualStagnation::ResidualStagnation( lama::NormPtr norm )
    : Criterion(), mNorm( norm ), mLookback( 1 ), //TODO check default value
      mLastResidualNorms( 1 ), mNextEntry( 0 ), mEntriesReady( false ), mPrecision( 0.1 ) //TODO check default value
{
}

ResidualStagnation::ResidualStagnation( lama::NormPtr norm, IndexType lookback, lama::Scalar precision )
    : Criterion(), mNorm( norm ), mLookback( lookback ), mLastResidualNorms( lookback ), mNextEntry( 0 ), mEntriesReady(
          false ), mPrecision( precision )
{
}

ResidualStagnation::ResidualStagnation( const ResidualStagnation& other )
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
    mLastResidualNorms[mNextEntry] = ( *mNorm )( solver.getResidual() );
    mNextEntry = ( mNextEntry + 1 ) % mLookback;

    if ( mNextEntry == 0 )
    {
        mEntriesReady = true;
    }

    if ( mEntriesReady )
    {
        lama::Scalar min = *std::min_element( mLastResidualNorms.begin(), mLastResidualNorms.end() );
        lama::Scalar max = *std::max_element( mLastResidualNorms.begin(), mLastResidualNorms.end() );
        min = std::max( std::numeric_limits<lama::Scalar>::min(), min );
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

const lama::NormPtr ResidualStagnation::getNorm() const
{
    return mNorm;
}

const lama::Scalar ResidualStagnation::getPrecision() const
{
    return mPrecision;
}

void ResidualStagnation::setLookback( IndexType lookback )
{
    mLookback = lookback;
}

void ResidualStagnation::setPrecision( const lama::Scalar precision )
{
    mPrecision = precision;
}

void ResidualStagnation::writeAt( std::ostream& stream ) const
{
    stream << "ResStgn<" << getPrecision() << ", " << getLookback() << ">";
}

} /* end namespace solver */

} /* end namespace scai */
