/**
 * @file AMGSetup.cpp
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
 * @brief AMGSetup.cpp
 * @author Jiri Kraus
 * @date 28.10.2011
 */

// hpp
#include <scai/solver/AMGSetup.hpp>

#include <scai/common/SCAITypes.hpp>
#include <scai/common/macros/instantiate.hpp>

namespace scai
{

namespace solver
{

SCAI_LOG_DEF_LOGGER( _AMGSetup::logger, "AMGSetup" )

/* ========================================================================= */
/*    static methods (for _AMGSetup - Factory, untyped )                     */
/* ========================================================================= */

_AMGSetup* _AMGSetup::getAMGSetup( const common::ScalarType scalarType, const std::string& setupType )
{
    return create( AMGSetupCreateKeyType( scalarType, setupType ) );
}

/* ========================================================================= */
/*    static methods (for AMGSetup<ValueType> - Factory                      */
/* ========================================================================= */

template<typename ValueType>
void AMGSetup<ValueType>::getCreateValues( std::vector<std::string>& values )
{
    std::vector<AMGSetupCreateKeyType> createValues;

    _AMGSetup::getCreateValues( createValues );  // all solvers ( valueType, solvertype )

    values.clear();

    for ( size_t i = 0; i < createValues.size(); ++i )
    {
        if ( createValues[i].first == common::TypeTraits<ValueType>::stype )
        {
            // AMGSetup for this value type
            values.push_back( createValues[i].second );
        }
    }
}

template<typename ValueType>
AMGSetup<ValueType>* AMGSetup<ValueType>::getAMGSetup( const std::string& setupType )
{
    _AMGSetup* setup = _AMGSetup::getAMGSetup( common::TypeTraits<ValueType>::stype, setupType );

    SCAI_ASSERT_DEBUG( dynamic_cast<AMGSetup<ValueType>*>( setup ), "Illegal setup" )

    return reinterpret_cast<AMGSetup<ValueType>*>( setup );
}

/* ========================================================================= */
/*    Constructor/Destructor                                                 */
/* ========================================================================= */

template<typename ValueType>
AMGSetup<ValueType>::AMGSetup() : 

    mHostOnlyLevel( std::numeric_limits<IndexType>::max() ), 
    mHostOnlyVars( 0 ), 
    mReplicatedLevel( std::numeric_limits<IndexType>::max() )
{
}

template<typename ValueType>
AMGSetup<ValueType>::~AMGSetup()
{
}

/* ========================================================================= */
/*    Methods                                                                */
/* ========================================================================= */

template<typename ValueType>
void AMGSetup<ValueType>::setHostOnlyLevel( IndexType hostOnlyLevel )
{
    mHostOnlyLevel = hostOnlyLevel;
}

template<typename ValueType>
void AMGSetup<ValueType>::setHostOnlyVars( IndexType hostOnlyVars )
{
    mHostOnlyLevel = hostOnlyVars;
}

template<typename ValueType>
void AMGSetup<ValueType>::setReplicatedLevel( IndexType replicatedLevel )
{
    mReplicatedLevel = replicatedLevel;
}

template<typename ValueType>
void AMGSetup<ValueType>::writeAt( std::ostream& stream ) const
{
    stream << "AMGSetup( ... )";
}

/* ========================================================================= */
/*       Template instantiations                                             */
/* ========================================================================= */

SCAI_COMMON_INST_CLASS( AMGSetup, SCAI_NUMERIC_TYPES_HOST )

} /* end namespace solver */

} /* end namespace scai */
