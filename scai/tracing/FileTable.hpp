/**
 * @file FileTable.hpp
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
