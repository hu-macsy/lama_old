/**
 * @file AMGSetupFactory.hpp
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
 * @brief Definition of a singleton class that provides an AMGSetup of registered
 *        AMGSetup types.
 *
 * @author Thomas Brandes
 * @date 15.03.2013
 * @since 1.0.0
 */

#ifndef LAMA_AMG_SETUP_FACTORY_HPP_
#define LAMA_AMG_SETUP_FACTORY_HPP_

// for dll_import
#include <common/config.hpp>

// base classes
#include <common/NonCopyable.hpp>

// others
#include <lama/solver/AMGSetup.hpp>

// logging
#include <logging.hpp>

// boost
#include <common/shared_ptr.hpp>

#include <map>

namespace lama
{

class AMGSetupManager;

/** @brief Singleton class for an object that provides AMG setups of a certain type.
 *
 *  The factory uses managers that must have been registered before. This gives
 *  extensibility for new AMG setup classes as the factory itself has no need
 *  to create any AMG setup by itself.
 *
 *  AMGSetup managers are stored via shared pointers so ownership of managers is
 *  flexible. The release function of the factory gives up ownership of all managers
 *  so they will be freed if there are no other references.
 *
 *  Note: The type of the AMG setup is given by a string.
 */

class COMMON_DLL_IMPORTEXPORT AMGSetupFactory: common::NonCopyable

{
public:

    typedef std::map<std::string,common::shared_ptr<AMGSetupManager> > AMGSetupToManagerMap;

    /** Get an AMGSetup of a certain type from the factory.
     *
     *  @param type is the name of the needed AMG setup.
     *  @returns pointer to the desired AMG setup, the default one if not found
     */

    static AMGSetupPtr get( const std::string& type );

    /** Get a default AMG setup from the factory.
     *
     *  @returns pointer to the default AMG setup.
     */

    static AMGSetupPtr get();

    /** Method that returns a reference to the AMG setup factory singleton. */

    static AMGSetupFactory& getFactory();

    /** Get a AMG setup of a certain type from the factory.
     *
     *  @param type is the name of the needed AMG setup.
     *  @returns pointer to the desired AMG setup, the default one if not found
     */
    common::shared_ptr<AMGSetupManager> getAMGSetupManager( const std::string& type );

    /**
     *  Releases all communication managers registered in the AMG setup factory instance.
     *  Might be called at the end of a program to finalize all communication before program exit.
     *  NOTE: was mandatory for VAMPIR trace of MPI communication.
     */
    static void release();

    /** This method adds a new AMG setup manager to the factory.
     *
     * @param type is the type of AMG setup that is managed by the manager.
     * @param manager is the AMG setup manager that provides AMGSetup of the given type.
     *
     * The factory takes ownership of the manager, i.e. the manager will be deleted if it is replaced or
     * if the factory is destroyed or released.
     */
    void addAMGSetupManager( const std::string& type, common::shared_ptr<AMGSetupManager> manager );

    /** Query routine for the default AMG setup type. */

    const std::string& getDefaultAMGSetupType() const;

    /** Set the default AMG setup type */

    void setDefaultAMGSetupType( const std::string& type ) const;

    /** Destructor, might also free all registered AMG setup managers. */

    virtual ~AMGSetupFactory();

private:

    /** Routine to find a default AMG setup. */

    void setDefaultAMGSetupType() const;

    /** @brief Constructor of a AMG setup factory.
     *
     *  The constructor is private to guarantee that only one singleton instance is created.
     */

    AMGSetupFactory();

    /** Map with manager for each registered AMG setup type. */

    AMGSetupToManagerMap mAMGSetupToManagerMap;

    mutable std::string mDefaultAMGSetupType; // name of the default AMG setup

    LAMA_LOG_DECL_STATIC_LOGGER( logger )
};

}

#endif // LAMA_AMG_SETUP_FACTORY_HPP_
