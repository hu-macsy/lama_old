/**
 * @file ContextAccess.hpp
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
 * @brief Definition of an access class to access a context.
 * @author Thomas Brandes
 * @date 14.07.2011
 */
#pragma once

// for dll_import
#include <scai/common/config.hpp>

// base classes
#include <scai/hmemo/Access.hpp>

// local library
#include <scai/hmemo/Context.hpp>

// logging
#include <scai/logging.hpp>

namespace scai
{

namespace hmemo
{

/** @brief This class is used for accessing a Context to enable computations on it.
 *
 *  This class guarantees that enable of a context is always followed by
 *  a disable of the context in case of exception.
 *
 *  Member variables contain source code location where access to the context
 *  was enabled.
 */
class COMMON_DLL_IMPORTEXPORT ContextAccess: public Access
{
public:

    /** Constructor for an access, enables the context. */

    ContextAccess( ContextPtr context, const char* filename, int line );

    /** Destructor for an access, disables the context. */

    ~ContextAccess();

    /** Mandatory routine for each access class. */

    virtual void release();

    virtual void writeAt( std::ostream& stream ) const;

private:

    SCAI_LOG_DECL_STATIC_LOGGER( logger )

    ContextAccess    (); // do not allow the default constructor

    const Context& mContext;//!<  context that is accessed.

    bool mReleased;//!< true, if access has already been released

    const char* mFile;//!<  file name of access

    int mLine;//!< line number of access
};

} /* end namespace hmemo */

} /* end namespace scai */

/** Macro that should be used for any context access
 *
 *  \code
 *  {
 *      SCAI_CONTEXT_ACCESS( cudaContext )
 *      ...   // computations on CUDA device
 *      // context is released at end of the scope by destructor of access object.
 *  }
 *  \endcode
 */

#define SCAI_CONTEXT_ACCESS( context ) scai::hmemo::ContextAccess SCAI_Ctx_Access__( context, __FILE__, __LINE__ );
