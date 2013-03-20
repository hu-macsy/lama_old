/**
 * @file AMGSetupManager.hpp
 *
 * @license
 * Copyright (c) 2011
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 * @endlicense
 *
 * @brief AMGSetupManager.hpp
 * @author Thomas Brandes
 * @date 15.03.2013
 * $Id$
 */

#ifndef LAMA_AMG_SETUP_MANAGER_HPP_
#define LAMA_AMG_SETUP_MANAGER_HPP_

// for dll_import
#include <lama/config.hpp>

// base classes
#include <lama/NonCopyable.hpp>

// others
#include <lama/solver/AMGSetup.hpp>

// logging
#include <logging/logging.hpp>

namespace lama
{

/** @brief Base class for a AMG setup manager.
 *
 *  For each AMG setup type a AMG setup manager must be available.
 *  The AMG setup manager registers in the AMGSetupFactory and
 *  delivers a shared pointer for a AMG setup.
 *
 */

class AMGSetupManager: private NonCopyable
{
public:

    virtual ~AMGSetupManager();

    /** @brief Method that returns a AMGSetup.  */

    virtual AMGSetupPtr getAMGSetup() = 0;

protected:

    /** Constructor of a AMGSetupManager, type must be specified. */

    AMGSetupManager( const char* type );

    std::string mAMGSetupType; //!< type of AMGSetup managed

    LAMA_LOG_DECL_STATIC_LOGGER(logger);
};

} // namespace lama

#endif // LAMA_AMG_SETUP_MANAGER_HPP_
