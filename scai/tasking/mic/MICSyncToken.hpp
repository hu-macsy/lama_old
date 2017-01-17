/**
 * @file MICSyncToken.hpp
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
 * @brief MICSyncToken.hpp
 * @author Thomas Brandes
 * @date 02.07.2013
 */

#pragma once

// base classes
#include <scai/tasking/SyncToken.hpp>

// internal scai libraries
#include <scai/common/shared_ptr.hpp>

namespace scai
{

namespace tasking
{

/** Class that sycnchronizes with a MIC offload transfer or computation. */

class COMMON_DLL_IMPORTEXPORT MICSyncToken: public SyncToken

{
public:

    /** Constructor for a MIC sychronization token.
     *
     *  @param[in]  device  is the MIC device where asynchronous operation takes place
     */

    MICSyncToken( int device );

    inline int getDevice() const
    {
        return mDevice;
    }

    inline int& signal()
    {
        return mSignal;
    }

    virtual ~MICSyncToken();

    /** After starting the offload computation/transfer with a signal this signal is set here. */

    void setSignal( int signal );

    virtual void wait();

    virtual bool probe() const;

    /** Get sync token in case of asynchronous execution should be started. */

    static MICSyncToken* getCurrentSyncToken();

private:

    int mSignal; // set by an offload computation
    int mDevice; // device set by constructor
};

} /* end namespace tasking */

} /* end namespace scai */
