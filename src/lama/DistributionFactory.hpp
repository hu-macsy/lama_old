/**
 * @file DistributionFactory.hpp
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
 * @brief Definition of a singleton class that provides distributions of registered
 *        distribution types.
 *
 * @author Thomas Brandes
 * @date 20.12.2012
 * $Id$
 */

#ifndef LAMA_DISTRIBUTION_FACTORY_HPP_
#define LAMA_DISTRIBUTION_FACTORY_HPP_

// for dll_import
#include <lama/config.hpp>

// base classes
#include <lama/NonCopyable.hpp>

// others
#include <lama/distribution/Distribution.hpp>

// logging
#include <logging/logging.hpp>

// boost
#include <boost/shared_ptr.hpp>

#include <map>

namespace lama
{

class DistributionManager;

/** @brief Singleton class for an object that provides distributions of a certain type.
 *
 *  The factory uses managers that must have been registered before. This gives
 *  extensibility for new distribution classes as the factory itself has no need
 *  to create any distribution by itself.
 *
 *  Distribution managers are stored via shared pointers so ownership of managers is
 *  flexible. The release function of the factory gives up ownership of all managers
 *  so they will be freed if there are no other references.
 *
 *  Note: The type of the disrtribution is given by a string.
 */

class LAMA_DLL_IMPORTEXPORT DistributionFactory: NonCopyable

{
public:

    typedef std::map<std::string,boost::shared_ptr<DistributionManager> > DistributionToManagerMap;

    /** @brief Get a distribution of a certain type with using additional arguments
     *
     *  Abbreviation for get( size, theFactory.getDefaultDistributionType(), theFactory.mDefaultArgs )
     */
    static DistributionPtr get( const IndexType size );

    /** @brief Get a distribution of a certain type with using additional arguments
     *
     *  Abbreviation for get( size, type, theFactory.mDefaultArgs )
     */
    static DistributionPtr get( const IndexType size, const std::string& type );

    /** @brief Get a distribution of a certain type with using additional arguments
     *
     *  Abbreviation for get( size, type, arguments )
     */
    static DistributionPtr get( const IndexType size, const std::string& type, int arg1 );

    /** @brief Get a distribution of a certain type with using additional arguments
     *
     *  Abbreviation for get( size, type, arguments )
     */
    static DistributionPtr get( const IndexType size, const std::string& type, int arg1, int arg2 );

    /** Get a distribution of a certain type with using additional arguments
     *
     * @param[in] size        is the global number of elements to be distributed
     * @param[in] type        is the communicator type needed
     * @param[in] arguments   additional arguments to the manager used for creation
     *
     * The additional arguments will be passed to the manager that uses them for creation
     */
    static DistributionPtr get( const IndexType size, const std::string& type, const std::vector<int>& arguments );

    /** Get a default distribution from the factory.
     *
     *  @returns pointer to the default distribution.
     */
    static DistributionFactory& getFactory();

    /** Method that sets the default distribution of the factory. */

    static void setDefaultDistribution( const std::string& type );
    static void setDefaultDistribution( const std::string& type, int arg1 );
    static void setDefaultDistribution( const std::string& type, int arg1, int arg2 );

    /** Get a communicator of a certain type from the factory.
     *
     *  @param type is the name of the needed communicator.
     *  @returns pointer to the desired communicator, the default one if not found
     */
    boost::shared_ptr<DistributionManager> getDistributionManager( const std::string& type );

    /**
     *  Releases all distribution managers registered in the communicator factory instance.
     *  Might be called at the end of a program to finalize all distribution before program exit.
     */
    static void release();

    /** This method adds a new distribution manager to the factory.
     *
     * @param type is the type of distribution that is managed by the manager.
     * @param manager is the distribution manager that provides Distribution of the given type.
     *
     * The factory takes ownership of the manager, i.e. the manager will be deleted if it is replaced or
     * if the factory is destroyed or released.
     */
    void addDistributionManager( const std::string& type, boost::shared_ptr<DistributionManager> manager );

    /** Query routine for the default distribution type. */

    const std::string& getDefaultDistributionType() const;

    /** Set the default distribution type */

    void setDefaultDistributionType( const std::string& type ) const;

    /** Destructor, also frees all registered communicator managers. */

    virtual ~DistributionFactory();

private:

    /** Routine to find a default communicator. */

    void setDefaultDistributionType() const;

    /** @brief Constructor of a communicator factory.
     *
     *  The constructor is private to guarantee that only one singleton instance is created.
     */

    DistributionFactory();

    /** Map with manager for each registered communicator type. */

    DistributionToManagerMap mDistributionToManagerMap;

    mutable std::string mDefaultDistributionType; // name of the default communicator

    mutable std::vector<int> mDefaultArgs;

    LAMA_LOG_DECL_STATIC_LOGGER( logger )
};

}

#endif // LAMA_DISTRIBUTION_FACTORY_HPP_
