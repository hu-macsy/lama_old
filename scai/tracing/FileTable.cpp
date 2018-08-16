/**
 * @file FileTable.cpp
 *
 * @license
 * Copyright (c) 2009-2018
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the SCAI framework LAMA.
 *
 * LAMA is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * LAMA is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
 * more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with LAMA. If not, see <http://www.gnu.org/licenses/>.
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
    MapFiles::iterator it = mapTimer.find( fileName );

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
