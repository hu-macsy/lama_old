/**
 * @file lama/io/IOUtils.cpp
 *
 * @license
 * Copyright (c) 2009-2016
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the Library of Accelerated Math Applications (LAMA).
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
 * @endlicense
 *
 * @brief ToDo: Missing description in ./lama/io/IOUtils.cpp
 * @author eschricker
 * @date 15.03.2016
 */

// base class
#include <scai/lama/io/IOUtils.hpp>

// scai internal libraries
#include <scai/common/macros/throw.hpp>

// others
#include <ostream>

namespace scai
{

namespace lama
{

IOUtils::file_size_t IOUtils::getFileSize( const char* filename )
{
    FILE* pFile = fopen( filename, "rb" );

    file_size_t size = 0;

    if( pFile == NULL )
    {
        COMMON_THROWEXCEPTION( "File " << filename << " could not be opened" )
    }
    else
    {
        fseek( pFile, 0, SEEK_END ); // non-portable
        size = ftell( pFile );
        fclose( pFile );
    }

    return size;
}

} /* end namespace lama */

} /* end namespace scai */
