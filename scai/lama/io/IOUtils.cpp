/*
 * IOUtils.cpp
 *
 *  Created on: Mar 15, 2016
 *      Author: eschricker
 */

// base class
#include <scai/lama/io/IOUtils.hpp>

// scai internal libraries
#include <scai/common/macros/throw.hpp>

// others
#include <ostream>

namespace scai {

namespace lama {

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
