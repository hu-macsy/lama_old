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
#include <scai/common/config.hpp>

// internal scai libraries
#include <scai/logging.hpp>

// std
#include <vector>
#include <cstring>
#include <string>
#include <map>
#include <cstdio>

namespace scai
{

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

    SCAI_LOG_DECL_STATIC_LOGGER( logger )

private:

    std::vector<std::string> array; //!<  Entries for all filenames

    typedef std::map<const std::string, int> MapFiles;

    /** Map of file strings to file ids that are the indexes to array.
     *
     *  std::string should be used as key and not const char* as this
     *  might cause problems in case of auto strings.
     */

    MapFiles mapTimer;
};

} /* end namespace tracing */

} /* end namespace scai */
