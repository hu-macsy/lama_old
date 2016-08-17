/**
 * @file CRTPCommunicator.hpp
 *
 * @license
 * Copyright (c) 2009-2016
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
 * @brief CRTP class that provides the polymorphism for the virtual routines and
 *        so derived classes have to provide only template routines
 * @author Thomas Brandes
 * @date 12.05.2014
 */

#pragma once

// for dll_import
#include <scai/common/config.hpp>

// base classes
#include <scai/dmemo/Communicator.hpp>

// internal scai libraris
#include <scai/tasking/NoSyncToken.hpp>

#include <scai/common/macros/loop.hpp>

// std
#include <vector>

namespace scai
{

namespace dmemo
{

/** This template class supports static polymorphism to define
 *  common routines for derived classes of Communicator
 *
 */

template<class Derived>
class COMMON_DLL_IMPORTEXPORT CRTPCommunicator: public Communicator
{
public:

#define SCAI_DMEMO_CRTP_COMMUNICATOR_METHODS( _type )                                       \
                                                                                            \
    virtual void maxloc( _type& val, IndexType& location, PartitionId root ) const          \
    {                                                                                       \
        static_cast<const Derived*>( this )->maxlocImpl( val, location, root );             \
    }                                                                                       \
                                                                                            \
    // define communicator methods for all supported data types

    SCAI_COMMON_LOOP( SCAI_DMEMO_CRTP_COMMUNICATOR_METHODS, SCAI_ALL_TYPES )

#undef SCAI_DMEMO_CRTP_COMMUNICATOR_METHODS

protected:

    // Default constructor can only be called by derived classes.

    CRTPCommunicator<Derived>( const CommunicatorKind& type )
        : Communicator( type )
    {
    }

};

} /* end namespace dmemo */

} /* end namespace scai */
