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
#include <iostream>

namespace scai
{

namespace benchmark
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

void InputSet::writeAt( std::ostream& stream ) const
{
    stream << "InputSet( id = " << mId << ", name = " << mName << " )";
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

InputSet* InputSet::parseAndCreate( const std::string& specification )
{
    std::string::size_type pos1 = specification.find_first_of( " (" );
    SCAI_ASSERT_NE_ERROR( pos1, std::string::npos, specification << ": no ( found" )
    std::string::size_type pos2 = specification.find_first_not_of( " (", pos1 );
    SCAI_ASSERT_NE_ERROR( pos2, std::string::npos, specification << ": no arg found" )

    std::string::size_type pos3 = specification.find_first_of( " )", pos2 );
    if ( pos3 == std::string::npos )
    {
        pos3 = specification.length();
    }

    std::string keyValue = specification.substr( 0, pos1 );
    std::string argument = specification.substr( pos2, pos3 - pos2  );

    std::cout << "pos1 = " << pos1 << ", pos2 = " << pos2 << ", pos3 = " << pos3 
              << ", keyValue=<" << keyValue << ">, arg=<" << argument << ">" << std::endl;

    return create( keyValue, argument );
}

}

}
