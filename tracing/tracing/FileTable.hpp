/**
 * @file FileTable.hpp
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
 * @brief Definition of class that maps file names to ids
 * @author Thomas Brandes
 * @date 11.06.2015
 */

#pragma once

// for dll_import
#include <common/config.hpp>

#include <logging/logging.hpp>

#include<vector>
#include<cstring>
#include<string>
#include<map>
#include<cstdio>

namespace tracing
{

/** Class that manages all used file names in a table */

class COMMON_DLL_IMPORTEXPORT FileTable
{

public:

    /** Constructor of a new file table.
     *
     *  The region table must contain the thread id as it might be 
     *  written later by main thread.
     */

    FileTable();

    /** Destructor. */

    ~FileTable();

    /** Get the id of a file, creates a new entry if file is not available yet.
     *
     *  @param[in] name   is the name of the region
     */

    int getFileId( const char* name );

    /** Get file name by its id. */

    const char* getFileName( int fileId );

protected:

    LAMA_LOG_DECL_STATIC_LOGGER( logger )

private:

    struct    CmpString
    {
        bool operator()( const char* a, const char* b ) const
        {
            return std::strcmp( a, b) < 0;
        }
    };

    std::vector<std::string> array; //!<  Entries for all filenames

    /** Map of file strings to file ids that are the indexes to array.
     *
     *  Timer strings are given by pointer; pointers will be always valid as string
     *  remains as member variable in array.
     */

    std::map<const char*, int, CmpString> mapTimer;
};

}
