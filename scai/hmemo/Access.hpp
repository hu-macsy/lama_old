/**
 * @file Access.hpp
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
 * @brief Definition of common base class for accesses to a HArray to enforce consistency.
 * @author Jiri Kraus
 * @date 20.05.2011
 */

#pragma once

// for dll_import
#include <scai/common/config.hpp>

// base classes
#include <scai/common/NonCopyable.hpp>
#include <scai/common/Printable.hpp>
#include <scai/tasking/SyncToken.hpp>

// logging
#include <scai/logging.hpp>

namespace scai
{

namespace hmemo
{

/**  class for all kind of access classes; accesses must always provide a release method. */

class COMMON_DLL_IMPORTEXPORT Access:

    public  common::Printable,
    private common::NonCopyable,
    public  tasking::SyncTokenMember
{
public:

    Access();

    virtual ~Access();

    virtual void release() = 0;

    virtual void writeAt( std::ostream& stream ) const;

protected:

    SCAI_LOG_DECL_STATIC_LOGGER( logger )
};

} /* end namespace hmemo */

} /* end namespace scai */
