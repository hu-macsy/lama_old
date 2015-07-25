/**
 * @file MICContextManager.hpp
 *
 * @license
 * Copyright (c) 2009-2015
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
 * @brief Definition of manager class to manage MIC contexts.
 * @author Thomas Brandes
 * @date 01.07.2013
 * @since 1.1.0
 */
#ifndef LAMA_MIC_CONTEXT_MANAGER_HPP_
#define LAMA_MIC_CONTEXT_MANAGER_HPP_

// for dll_import
#include <common/config.hpp>

// base classes
#include <lama/ContextManager.hpp>

// others
#include <lama/Context.hpp>

// boost
#include <common/weak_ptr.hpp>

#include <string>

namespace lama
{

class MICContext;

/** maximal number of MIC devices that might be available */

#define LAMA_MAX_MIC_DEVICES 4
#define LAMA_MIC_ENV_FOR_DEVICE "LAMA_MIC_DEVICE"

/** @brief This class is a derived class of ContextManager for managing MIC Context.
 *
 *  It provides instances of MICContext for a given device number.
 */

class COMMON_DLL_IMPORTEXPORT MICContextManager: public ContextManager
{
public:

    /** Implementation of ContextManager::getContext */

    ContextPtr getContext( int deviceNr )
    {
        return getInstance( deviceNr );
    }

    /** Destructor will free all managed contexts. */

    ~MICContextManager();

    /** Static version of the member method getContext. */

    static ContextPtr getInstance( int deviceNr );

private:

    /** Make constructor private, only single static object will be created. */

    MICContextManager();

    /** The manager keeps weak pointers to MIC contexts so it can return
     *  a shared pointer if the contexts are still in use.
     *  Entry is available for each possible device.
     */

    static boost::weak_ptr<MICContext> mMICContext[LAMA_MAX_MIC_DEVICES];

    static MICContextManager theInstance; //!< singleton instance of this class

    static int defaultDeviceNr; //!< will be set to the default device

    /** Getter for the default device number; it will be set if not available at first call. */

    static int getDefaultDeviceNr(); // getter for defaultDeviceNr, might be set
};

}

#endif // LAMA_MIC_CONTEXT_MANAGER_HPP_
