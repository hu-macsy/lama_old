/**
 * @file CollectiveFile.cpp
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
 * @brief Collectice/concurrent file I/O 
 * @author Thomas Brandes
 * @date 14.01.2019
 */

#include <scai/dmemo/CollectiveFile.hpp>

namespace scai
{

namespace dmemo
{

/* -------------------------------------------------------------------------- */

SCAI_LOG_DEF_LOGGER( CollectiveFile::logger, "CollectiveFile" )

/* -------------------------------------------------------------------------- */

CollectiveFile::CollectiveFile( CommunicatorPtr comm ) :

    mComm( comm ),
    mOffset( 0 )
{
}

CollectiveFile::~CollectiveFile()
{
}

/* -------------------------------------------------------------------------- */

void CollectiveFile::set( const char* filename, size_t offset )
{
    mFileName = filename;
    mOffset = offset;
}

/* -------------------------------------------------------------------------- */

}

}
