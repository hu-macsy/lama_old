/**
 * @file ResidualThreshold.cpp
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

#include <scai/common/macros/instantiate.hpp>

namespace scai
{

namespace solver
{

std::ostream& operator<<( std::ostream& stream, const ResidualCheck checkMode )
{
    switch ( checkMode )
    {
        case ResidualCheck::Relative:
            stream << "Realitve";
            break;
        case ResidualCheck::Absolute:
            stream << "Absolute";
            break;
        case ResidualCheck::Divergence:
            stream << "Divergence";
            break;
        default:
            stream << "<unknown_residual_threshold_check_mode";
    }

    return stream;
}

template<typename ValueType>
ResidualThreshold<ValueType>::ResidualThreshold() : 

    Criterion<ValueType>(), 
    mNorm( new lama::L2Norm<ValueType>() ),
    mCheckMode( ResidualCheck::Relative ), 
    mPrecision( ValueType( 1e-5 ) ), 
    mFirstNormResult( ValueType( -1 ) )
{
}

template<typename ValueType>
ResidualThreshold<ValueType>::ResidualThreshold( const lama::NormPtr<ValueType> norm ) : 

    Criterion<ValueType>(), 
    mNorm( norm ), 
    mCheckMode( ResidualCheck::Relative ), 
    mPrecision( ValueType( 1e-5 ) ), 
    mFirstNormResult( ValueType( -1 ) )
{
    SCAI_ASSERT_ERROR( mNorm, "norm must not be NULL pointer" )
}

template<typename ValueType>
ResidualThreshold<ValueType>::ResidualThreshold( const lama::NormPtr<ValueType> norm, ValueType precision, ResidualCheck checkMode ) : 

    Criterion<ValueType>(), 
    mNorm( norm ), 
    mCheckMode( checkMode ), 
    mPrecision( precision ), 
    mFirstNormResult( -1.0 )
{
    SCAI_ASSERT_ERROR( mNorm, "norm must not be NULL pointer" )
}

template<typename ValueType>
ResidualThreshold<ValueType>::ResidualThreshold( const ResidualThreshold<ValueType>& other ) : 

    Criterion<ValueType>(), 
    mNorm( other.mNorm ), 
    mCheckMode( other.mCheckMode ), 
    mPrecision( other.mPrecision ), 
    mFirstNormResult( other.mFirstNormResult )
{
}

template<typename ValueType>
ResidualThreshold<ValueType>::~ResidualThreshold()
{
}

template<typename ValueType>
ResidualThreshold<ValueType>* ResidualThreshold<ValueType>::copy() const
{
    return new ResidualThreshold<ValueType>( *this );
}

template<typename ValueType>
inline bool ResidualThreshold<ValueType>::isSatisfied( const IterativeSolver<ValueType>& solver )
{
    typedef RealType<ValueType> DefaultReal;

    DefaultReal normResult = ( *mNorm )( solver.getResidual() );
    SCAI_ASSERT( normResult >= DefaultReal( 0 ), "A norm should be always positive but is " << normResult );

    switch ( mCheckMode )
    {
        case ResidualCheck::Absolute:
            SCAI_LOG_DEBUG( logger,
                            "Absolute residual in iteration " << solver.getIterationCount() 
                            << " is " << normResult << " should become smaller than precision " << mPrecision )

            return normResult < mPrecision;

        case ResidualCheck::Relative:
        {
            if ( mFirstNormResult == ValueType( -1 ) )
            {
                //TODO define member variable for solver with getInitialResidual function
                mFirstNormResult = normResult;
            }

            SCAI_LOG_DEBUG( logger,
                            "Relative residual in iteration " << solver.getIterationCount() 
                             << " is " << normResult << " divided by firstNormResult " << mFirstNormResult 
                             << " is " << normResult / mFirstNormResult 
                             << " should become smaller than precision " << mPrecision );

            return ( normResult / mFirstNormResult ) < mPrecision;
        }

        case ResidualCheck::Divergence:
        {
            if ( mFirstNormResult == -1.0 )
            {
                //TODO define member variable for solver with getInitialResidual function
                mFirstNormResult = normResult;
            }

            SCAI_LOG_DEBUG( logger,
                            "Relative residual in iteration " << solver.getIterationCount() << " is " << normResult << " divided by firstNormResult " << mFirstNormResult << " is " << normResult / mFirstNormResult << " should become larger than precision " << mPrecision );
            return ( normResult / mFirstNormResult ) > mPrecision;
        }
    }

    return false;
}

template<typename ValueType>
ResidualCheck ResidualThreshold<ValueType>::getCheckMode() const
{
    return mCheckMode;
}

template<typename ValueType>
ValueType ResidualThreshold<ValueType>::getFirstNormResult() const
{
    return mFirstNormResult;
}

template<typename ValueType>
const lama::NormPtr<ValueType> ResidualThreshold<ValueType>::getNorm() const
{
    return mNorm;
}

template<typename ValueType>
ValueType ResidualThreshold<ValueType>::getPrecision() const
{
    return mPrecision;
}

template<typename ValueType>
void ResidualThreshold<ValueType>::setCheckMode( ResidualCheck checkMode )
{
    mCheckMode = checkMode;
}

template<typename ValueType>
void ResidualThreshold<ValueType>::setFirstNormResult( ValueType firstNormResult )
{
    mFirstNormResult = firstNormResult;
}

template<typename ValueType>
void ResidualThreshold<ValueType>::setPrecision( ValueType precision )
{
    mPrecision = precision;
}

template<typename ValueType>
void ResidualThreshold<ValueType>::writeAt( std::ostream& stream ) const
{
    stream << "ResidualThreshold<" << getPrecision() << ", " << mCheckMode << ">";
}

/* ========================================================================= */
/*       Template instantiations                                             */
/* ========================================================================= */

SCAI_COMMON_INST_CLASS( ResidualThreshold, SCAI_NUMERIC_TYPES_HOST )

} /* end namespace solver */

} /* end namespace scai */
