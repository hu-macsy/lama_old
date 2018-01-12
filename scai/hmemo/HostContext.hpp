/**
 * @file HostContext.hpp
 *
 * @license
 * Copyright (c) 2009-2017
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the SCAI framework LAMA.
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
 *
 * Other Usage
 * Alternatively, this file may be used in accordance with the terms and
 * conditions contained in a signed written agreement between you and
 * Fraunhofer SCAI. Please contact our distributor via info[at]scapos.com.
 * @endlicense
 *
 * @brief Definition of context class that stands for execution on the host,
 *        i.e. convenctional CPU or multiple CPUs.
 * @author Thomas Brandes
 * @date 10.07.2011
 */

#pragma once

// for dll_import
#include <scai/common/config.hpp>

// base classes
#include <scai/hmemo/Context.hpp>

// internal scai libraries
#include <scai/tasking/TaskSyncToken.hpp>

#include <memory>

namespace scai
{

namespace hmemo
{

/** @brief This class implements the default HOST context.
 *
 *  This class is implemented as a singleton, only one default host
 *  context is available.
 *
 *  The default host context allocates/frees data in the usual way.
 */

class COMMON_DLL_IMPORTEXPORT HostContext:

    public Context,
    private Context::Register<HostContext>,
    public std::enable_shared_from_this<HostContext>
{

public:

    virtual ~HostContext();

    /** Static method required for create to use in Context::Register */

    static ContextPtr create( int deviceNr );

    /** Static method required for Context::Register */

    static common::ContextType createValue();

    /** Override Printable::writeAt with version for this class. */

    virtual void writeAt( std::ostream& stream ) const;

    /** Implement pure method Context::canUseMemory */

    virtual bool canUseMemory( const Memory& memory ) const;

    virtual tasking::TaskSyncToken* getSyncToken() const;

    virtual MemoryPtr getLocalMemoryPtr() const;

private:

    HostContext();

    SCAI_LOG_DECL_STATIC_LOGGER( logger )
};

inline common::ContextType HostContext::createValue()
{
    return common::ContextType::Host;
}

} /* end namespace hmemo */

} /* end namespace scai */
