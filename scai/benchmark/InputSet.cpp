/**
 * @file InputSet.cpp
 *
 * @license
 * Copyright (c) 2011
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
 * @brief InputSet.cpp
 * @author jiri
 * @date 06.04.2011
 * $Id$
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
