/**
 * @file common/thread.hpp
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
 * @brief Definition of functions to name threads.
 * @author Thomas Brandes
 * @date 11.06.2015
 */
#pragma once

// for dll_import
#include <scai/common/config.hpp>

// local library

// C++11 features

#include <memory>
#include <thread>
#include <mutex>
#include <condition_variable>

namespace scai
{

namespace common
{

/** 
 * Namespace that provides some free functions to give threads a name.
 *
 * Note: query a name returns a shared pointer to a string to avoid race conditions
 */
namespace thread
{
    typedef std::thread::id Id;

    /** Set a name for the current thread. */

    void defineCurrentThreadName( const char* name );

    /** Query the name of a thread. */

    std::shared_ptr<std::string> getThreadName( Id id );

    /** Query the name of the current thread. */

    std::shared_ptr<std::string> getCurrentThreadName();

} /* end namesapce thread */

} /* end namespace common */

} /* end namespace scai */

