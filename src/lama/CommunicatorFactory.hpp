/**
 * @file CommunicatorFactory.hpp
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
 * @brief Definition of a singleton class that provides communicators of registered
 *        communicator types.
 *
 * @author Jiri Kraus, Thomas Brandes
 * @date 23.02.2011
 * @since 1.0.0
 */

#ifndef LAMA_COMMUNICATOR_FACTORY_HPP_
#define LAMA_COMMUNICATOR_FACTORY_HPP_

// for dll_import
#include <lama/config.hpp>

// base classes
#include <lama/NonCopyable.hpp>

// others
#include <lama/Communicator.hpp>

// logging
#include <logging/logging.hpp>

// boost
#include <boost/shared_ptr.hpp>

#include <map>

namespace lama
{

class CommunicatorManager;

/** @brief Singleton class for an object that provides communicators of a certain type.
 *
 *  The factory uses managers that must have been registered before. This gives
 *  extensibility for new communicator classes as the factory itself has no need
 *  to create any communicator by itself.
 *
 *  Communicator managers are stored via shared pointers so ownership of managers is
 *  flexible. The release function of the factory gives up ownership of all managers
 *  so they will be freed if there are no other references.
 *
 *  Note: The type of the communicator is given by a string.
 */

class LAMA_DLL_IMPORTEXPORT CommunicatorFactory: NonCopyable

{
public:

    typedef std::map<std::string,boost::shared_ptr<CommunicatorManager> > CommunicatorToManagerMap;

    /** Get a communicator of a certain type with using inline arguments.
     *
     * @param type is the communicator type needed
     * @param argc pointer to the number of command line arguments
     * @param argv reference to the array of command line arguments
     *
     * Note: the command line arguments might be modified.
     */
    static CommunicatorPtr get( const std::string& type, int& argc, char** & argv );

    /** Get a communicator of a certain type from the factory.
     *
     *  @param type is the name of the needed communicator.
     *  @returns pointer to the desired communicator, the default one if not found
     */

    static CommunicatorPtr get( const std::string& type );

    /** Get a default communicator from the factory.
     *
     *  @returns pointer to the default communicator.
     */

    static CommunicatorPtr get();

    /** Method that returns a reference to the communicator factory singleton. */

    static CommunicatorFactory& getFactory();

    /** Get a communicator of a certain type from the factory.
     *
     *  @param type is the name of the needed communicator.
     *  @returns pointer to the desired communicator, the default one if not found
     */
    boost::shared_ptr<CommunicatorManager> getCommunicatorManager( const std::string& type );

    /**
     *  Releases all communication managers registered in the communicator factory instance.
     *  Might be called at the end of a program to finalize all communication before program exit.
     *  NOTE: was mandatory for VAMPIR trace of MPI communication.
     */
    static void release();

    /** This method adds a new communicator manager to the factory.
     *
     * @param type is the type of communicator that is managed by the manager.
     * @param manager is the communicator manager that provides Communicator of the given type.
     *
     * The factory takes ownership of the manager, i.e. the manager will be deleted if it is replaced or
     * if the factory is destroyed or released.
     */
    void addCommunicatorManager( const std::string& type, boost::shared_ptr<CommunicatorManager> manager );

    /** Query routine for the default communicator type. */

    const std::string& getDefaultCommunicatorType() const;

    /** Set the default communicator type */

    void setDefaultCommunicatorType( const std::string& type ) const;

    /** Destructor, also frees all registered communicator managers. */

    virtual ~CommunicatorFactory();

private:

    /** Routine to find a default communicator. */

    void setDefaultCommunicatorType() const;

    /** @brief Constructor of a communicator factory.
     *
     *  The constructor is private to guarantee that only one singleton instance is created.
     */

    CommunicatorFactory();

    /** Map with manager for each registered communicator type. */

    CommunicatorToManagerMap mCommunicatorToManagerMap;

    mutable std::string mDefaultCommunicatorType; // name of the default communicator

    LAMA_LOG_DECL_STATIC_LOGGER( logger )};

}

#endif // LAMA_COMMUNICATOR_FACTORY_HPP_
