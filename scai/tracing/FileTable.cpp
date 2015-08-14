/**
 * @file FileTable.cpp
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
 * @brief Implementation of methods for class FileTable.
 * @author Thomas Brandes
 * @date 11.06.2015
 */

// hpp
#include <scai/tracing/FileTable.hpp>

using namespace std;

namespace scai
{

namespace tracing
{

/* ---------------------------------------------------------------------- */

SCAI_LOG_DEF_LOGGER( FileTable::logger, "FileTable" )

/* ---------------------------------------------------------------------- */

FileTable::FileTable()
{
    SCAI_LOG_DEBUG( logger, "Constructor FileTable" )
    // avoid too much reallocations at the beginning
    array.reserve( 16 );
}

/* ---------------------------------------------------------------------- */

FileTable::~FileTable()
{
    SCAI_LOG_DEBUG( logger, "~FileTable" )
}

/* ---------------------------------------------------------------------- */

int FileTable::getFileId( const char* fileName )
{
    std::map<const char*, int, CmpString>::iterator it = mapTimer.find( fileName );

    if ( it == mapTimer.end() )
    {
        int fileId = array.size();
        array.push_back( fileName );
        const std::string& name = array[fileId];
        mapTimer.insert( std::pair<const char*, int>( name.c_str(), fileId ) );
        return fileId;
    }
    else
    {
        // alread defined
        return it->second;
    }
}

/* ---------------------------------------------------------------------- */

const char* FileTable::getFileName( int fileId )
{
    return array[fileId].c_str();
}

} /* end namespace tracing */

} /* end namespace scai */
