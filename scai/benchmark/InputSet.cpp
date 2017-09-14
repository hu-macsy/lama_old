/**
 * @file InputSet.cpp
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
 * @brief InputSet.cpp
 * @author jiri
 * @date 06.04.2011
 */
/**
 * @file InputSet.cpp
 * @author jiri
 * Created on: 04.05.2010
 */

#include <scai/benchmark/InputSet.hpp>
#include <utility>

namespace bf
{

InputSet::InputSet( const std::string& id )
    : mId( id ), mName( id )
{
}

InputSet::InputSet( const std::string& id, const std::string& name )
    : mId( id ), mName( name )
{
}

InputSet::~InputSet()
{
}

const std::string& InputSet::getId() const
{
    return mId;
}

const std::string& InputSet::getName() const
{
    return mName;
}

unsigned long InputSet::getNumFloatingPointOperations( const std::string& gid ) const
{
    ComplexityMap::const_iterator elem = mFlopMap.find( gid );
    if( elem == mFlopMap.end() )
    {
        return 0;
    }
    return elem->second;
}

void InputSet::setNumFloatingPointOperations( const std::string& gid, const unsigned long numFlops )
{
    mFlopMap[gid] = numFlops;
}

unsigned long InputSet::getProcessedBytes( const std::string& gid, const unsigned short sizeOfValueType ) const
{
    std::map<unsigned short,ComplexityMap>::const_iterator mapElem = mBWMap.find( sizeOfValueType );
    if( mapElem == mBWMap.end() )
    {
        return 0;
    }
    ComplexityMap::const_iterator elem = mapElem->second.find( gid );
    if( elem == mapElem->second.end() )
    {
        return 0;
    }
    return elem->second;
}

void InputSet::setProcessedBytes(
    const std::string& gid,
    const unsigned short sizeOfValueType,
    const unsigned long numBytes )
{
    mBWMap[sizeOfValueType][gid] = numBytes;
}

}
