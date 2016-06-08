/**
 * @file TaskSyncToken.hpp
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
 * @endlicense
 *
 * @brief Definition of a class that starts new Task.
 * @author Thomas Brandes
 * @date 04.05.2011
 */

#pragma once

// for dll_import
#include <scai/common/config.hpp>

// base classes
#include <scai/tasking/SyncToken.hpp>

// local library
#include <scai/tasking/Task.hpp>

// internal scai libraries
#include <scai/common/shared_ptr.hpp>

namespace scai
{

namespace tasking
{

/** Class to run a function asynchronously as a asynchronous Task. */

class COMMON_DLL_IMPORTEXPORT TaskSyncToken: public SyncToken
{
public:

    /** This constructor starts a routine as a Task and
     *  gives the SyncToken to wait on the completion.
     */
    TaskSyncToken( common::function<void()> routine, int numOmpThreads = 0 );

    /** This constructor starts a function as a Task and
     *  gives the SyncToken to wait on the completion.
     */
    TaskSyncToken();

    virtual void run( common::function<void()> function, int numOmpThreads = 0 );

    virtual ~TaskSyncToken();

    /** Wait for the completion of the Task. */

    virtual void wait();

    /** Test for the completion of the Task. */

    virtual bool probe() const;

    virtual void writeAt( std::ostream& stream ) const;

    /** Get sync token in case of asynchronous execution should be started. */

    static TaskSyncToken* getCurrentSyncToken();

private:

    common::shared_ptr<tasking::Task> mTask;

    SCAI_LOG_DECL_STATIC_LOGGER( logger )
};

} /* end namespace tasking */

} /* end namespace scai */
