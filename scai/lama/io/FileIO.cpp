/**
 * @file FileIO.cpp
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
 * @brief Implementation of static IO routines for matrix storage
 * @author Thomas Brandes
 * @date 10.06.2013
 * @since 1.0.1
 */

// hpp
#include <scai/lama/io/FileIO.hpp>

#include <scai/common/Assert.hpp>


namespace scai
{

namespace lama
{

// Help function to get the size of a file

FileIO::file_size_t FileIO::getFileSize( const char* filename )
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
