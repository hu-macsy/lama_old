/**
 * @file common/Thread.hpp
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
 * @brief Definition of the class Thread to manage thread ids and names
 * @author Thomas Brandes
 * @date 11.06.2015
 */
#pragma once

// for dll_import
#include <scai/common/config.hpp>

// local library
#include <scai/common/NonCopyable.hpp>
#include <scai/common/macros/throw.hpp>

// C++11 features

#include <memory>
#include <thread>
#include <mutex>
#include <condition_variable>

namespace scai
{

namespace common
{

/** This class helps to deal with thread ids and to name threads.
 *  It also supports critical regions by mutex and scoped locks.
 */

class COMMON_DLL_IMPORTEXPORT Thread : common::NonCopyable
{
public:

    typedef std::thread::id Id;

    /** Set a name for the current thread. */

    static void defineCurrentThreadName( const char* name );

    /** Query the name of a thread. */

    static const char* getThreadName( Id id );

    /** Query the name of the current thread. */

    static const char* getCurrentThreadName();

private:

    Thread();

};

} /* end namespace common */

} /* end namespace scai */

