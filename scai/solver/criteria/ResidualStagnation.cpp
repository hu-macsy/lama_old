/**
 * @file solver/criteria/ResidualStagnation.cpp
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
 * @brief ResidualStagnation.hpp
 * @author Kai Buschulte
 * @date 25.07.2011
 */

// hpp
#include <scai/solver/criteria/ResidualStagnation.hpp>

// local library
#include <scai/solver/IterativeSolver.hpp>
#include <scai/common/macros/instantiate.hpp>

// std
#include <iostream>
#include <algorithm>

namespace scai
{

namespace solver
{

template<typename ValueType>
ResidualStagnation<ValueType>::ResidualStagnation( lama::NormPtr<ValueType> norm ) : 

    Criterion<ValueType>(), 
    mNorm( norm ), 
    mLookback( 1 ), //TODO check default value
    mLastResidualNorms( 1 ), 
    mNextEntry( 0 ), 
    mEntriesReady( false ), 
    mPrecision( 0.1 ) //TODO check default value
{
}

template<typename ValueType>
ResidualStagnation<ValueType>::ResidualStagnation( lama::NormPtr<ValueType> norm, IndexType lookback, ValueType precision ) : 

    Criterion<ValueType>(), 
    mNorm( norm ), 
    mLookback( lookback ), 
    mLastResidualNorms( lookback ), 
    mNextEntry( 0 ), 
    mEntriesReady( false ), 
    mPrecision( precision )
{
}

template<typename ValueType>
ResidualStagnation<ValueType>::ResidualStagnation( const ResidualStagnation& other ) : 

    Criterion<ValueType>(), 
    mNorm( other.mNorm ), 
    mLookback( other.mLookback ), 
    mLastResidualNorms( other.mLastResidualNorms ), 
    mNextEntry( other.mNextEntry ), 
    mEntriesReady( other.mEntriesReady ), 
    mPrecision( other.mPrecision )
{
}

template<typename ValueType>
ResidualStagnation<ValueType>::~ResidualStagnation()
{
}

template<typename ValueType>
ResidualStagnation<ValueType>* ResidualStagnation<ValueType>::copy() const
{
    return new ResidualStagnation<ValueType>( *this );
}

template<typename ValueType>
bool ResidualStagnation<ValueType>::isSatisfied( const IterativeSolver<ValueType>& solver )
{
    mLastResidualNorms[mNextEntry] = ( *mNorm )( solver.getResidual() );
    mNextEntry = ( mNextEntry + 1 ) % mLookback;

    if ( mNextEntry == 0 )
    {
        mEntriesReady = true;
    }

    if ( mEntriesReady )
    {
        RealType<ValueType> min = *std::min_element( mLastResidualNorms.begin(), mLastResidualNorms.end() );
        RealType<ValueType> max = *std::max_element( mLastResidualNorms.begin(), mLastResidualNorms.end() );
        min = std::max( std::numeric_limits<RealType<ValueType> >::min(), min );
        //std::cout<< " max ="<<max<<"       min = "<<min<<"      max/min = "<<max/min<<"1+p = "<<(1.0+mPrecision)<<std::endl;
        mEntriesReady = false;
        return ( ( max / min ) < ( RealType<ValueType>( 1 ) + mPrecision ) );
    }

    return false;
}

template<typename ValueType>
IndexType ResidualStagnation<ValueType>::getLookback() const
{
    return mLookback;
}

template<typename ValueType>
const lama::NormPtr<ValueType> ResidualStagnation<ValueType>::getNorm() const
{
    return mNorm;
}

template<typename ValueType>
const ValueType ResidualStagnation<ValueType>::getPrecision() const
{
    return mPrecision;
}

template<typename ValueType>
void ResidualStagnation<ValueType>::setLookback( IndexType lookback )
{
    mLookback = lookback;
}

template<typename ValueType>
void ResidualStagnation<ValueType>::setPrecision( const ValueType precision )
{
    mPrecision = precision;
}

template<typename ValueType>
void ResidualStagnation<ValueType>::writeAt( std::ostream& stream ) const
{
    stream << "ResidualStagnation<" << getPrecision() << ", " << getLookback() << ">";
}

/* ========================================================================= */
/*       Template instantiations                                             */
/* ========================================================================= */

SCAI_COMMON_INST_CLASS( ResidualStagnation, SCAI_NUMERIC_TYPES_HOST )

} /* end namespace solver */

} /* end namespace scai */
