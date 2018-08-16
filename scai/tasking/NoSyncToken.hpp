/**
 * @file NoSyncToken.hpp
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
 * @brief Class that provides a sync token for complete synchronous execution.
 * @author Thomas Brandes
 * @date 25.03.2011
 */

#pragma once

// for dll_import
#include <scai/common/config.hpp>

// base classes
#include <scai/tasking/SyncToken.hpp>

namespace scai
{

namespace tasking
{

class COMMON_DLL_IMPORTEXPORT NoSyncToken: public SyncToken
{
public:

    inline NoSyncToken();

    inline virtual ~NoSyncToken();

    inline virtual void wait();

    inline virtual bool probe() const;

    inline void writeAt( std::ostream& stream ) const;
};

NoSyncToken::NoSyncToken()
{
    SCAI_LOG_DEBUG( logger, "NoSyncToken constructed" )
}

NoSyncToken::~NoSyncToken()
{
    SCAI_LOG_DEBUG( logger, "~NoSyncToken, synchronized = " << isSynchronized() )

    if ( !isSynchronized() )
    {
        setSynchronized(); // Important: accesses should be freed
    }
}

void NoSyncToken::writeAt( std::ostream& stream ) const
{
    stream << "NoSyncToken( synchronized = " << isSynchronized() << " )";
}

void NoSyncToken::wait()
{
    if ( !isSynchronized() )
    {
        setSynchronized(); // Important: accesses should be freed
    }
}

bool NoSyncToken::probe() const
{
    return true; // always ready
}

} /* end namespace tasking */

} /* end namespace scai */
