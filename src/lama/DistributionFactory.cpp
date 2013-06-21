/**
 * @file DistributionFactory.cpp
 *
 * @license
 * Copyright (c) 2009-2013
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
 * @brief DistributionFactory.cpp
 * @author Thomas Brandes, Jiri Kraus
 * @date 23.02.2011
 * @since 1.0.0
 */

// hpp
#include <lama/DistributionFactory.hpp>

// others
#include <lama/DistributionManager.hpp>
#include <lama/exception/Exception.hpp>
#include <lama/Settings.hpp>

// logging
#include <logging/logging.hpp>

// boost
#include <boost/shared_ptr.hpp>

using namespace std;
using namespace boost;

/** Name of the environment variable for the default distribution. */

#define LAMA_ENV_FOR_DISTRIBUTION "LAMA_DISTRIBUTION"

namespace lama
{

class Distribution;

/* ---- Static variables of DistributionFactory  ----------------------------------- */

LAMA_LOG_DEF_LOGGER( DistributionFactory::logger, "DistributionFactory" )

DistributionFactory& DistributionFactory::getFactory()
{
    // construct on first use idiom

    static std::auto_ptr<DistributionFactory> theDistributionFactory;

    // Note: using the auto_ptr guarantees that the instance will be freed at end of the program

    if ( !theDistributionFactory.get() )
    {
        // printf("create the communicator factory instance\n");

        theDistributionFactory = std::auto_ptr<DistributionFactory>( new DistributionFactory() );
    }

    return *theDistributionFactory;
}

DistributionPtr DistributionFactory::get(
    const IndexType size,
    const std::string& type,
    const std::vector<int>& arguments )
{
    LAMA_LOG_DEBUG( logger, "Get distribution of type " << type )

    const DistributionFactory& factory = getFactory();

    DistributionToManagerMap::const_iterator id = factory.mDistributionToManagerMap.find( type );

    if ( id != factory.mDistributionToManagerMap.end() )
    {
        return id->second->getDistribution( size, arguments );
    }

    const std::string& defaultType = factory.getDefaultDistributionType();

    LAMA_LOG_WARN( logger, "Distribution for " << type << " not supported, will use " << defaultType )

    // try to get a default communicator

    id = factory.mDistributionToManagerMap.find( factory.getDefaultDistributionType() );

    if ( id != factory.mDistributionToManagerMap.end() )
    {
        return id->second->getDistribution( size, arguments );
    }

    LAMA_THROWEXCEPTION( "Distribution for " << type << " not supported, no manager registered" )
}

DistributionPtr DistributionFactory::get( const IndexType size )
{
    const DistributionFactory& theFactory = getFactory();

    return get( size, theFactory.getDefaultDistributionType(), theFactory.mDefaultArgs );
}

DistributionPtr DistributionFactory::get( const IndexType size, const std::string& type )
{
    std::vector<int> arguments;

    return get( size, type, arguments );
}

DistributionPtr DistributionFactory::get( const IndexType size, const std::string& type, int arg1 )
{
    std::vector<int> arguments;
    arguments.push_back( arg1 );

    return get( size, type, arguments );
}

DistributionPtr DistributionFactory::get( const IndexType size, const std::string& type, int arg1, int arg2 )
{
    std::vector<int> arguments;
    arguments.push_back( arg1 );
    arguments.push_back( arg2 );

    return get( size, type, arguments );
}

/* ---- set default distribution --------------------------------------------------- */

void DistributionFactory::setDefaultDistribution( const std::string& type )
{
    const DistributionFactory& theFactory = getFactory();

    theFactory.setDefaultDistributionType( type );

    theFactory.mDefaultArgs.clear();
}

void DistributionFactory::setDefaultDistribution( const std::string& type, const int arg1 )
{
    const DistributionFactory& theFactory = getFactory();

    theFactory.setDefaultDistributionType( type );

    theFactory.mDefaultArgs.clear();
    theFactory.mDefaultArgs.push_back( arg1 );
}

void DistributionFactory::setDefaultDistribution( const std::string& type, const int arg1, const int arg2 )
{
    const DistributionFactory& theFactory = getFactory();

    theFactory.setDefaultDistributionType( type );

    theFactory.mDefaultArgs.clear();
    theFactory.mDefaultArgs.push_back( arg1 );
    theFactory.mDefaultArgs.push_back( arg2 );
}

const std::string& DistributionFactory::getDefaultDistributionType() const
{
    if ( mDefaultDistributionType.length() == 0 )
    {
        // default distribution has not been defined yet, so do it

        std::string defaultDistribution = "NO";

        // defaultDistribution might be overwritten by environment variable

        Settings::getEnvironment( defaultDistribution, LAMA_ENV_FOR_DISTRIBUTION );

        setDefaultDistributionType( defaultDistribution );
    }

    return mDefaultDistributionType;
}

boost::shared_ptr<DistributionManager> DistributionFactory::getDistributionManager( const std::string& type )
{
    LAMA_ASSERT( mDistributionToManagerMap.find( type ) != mDistributionToManagerMap.end(),
                 "Distribution for " << type << " not supported, no manager registered" )
    return mDistributionToManagerMap[type];
}

void DistributionFactory::setDefaultDistributionType( const std::string& defaultType ) const
{
    if ( mDistributionToManagerMap.find( defaultType ) != mDistributionToManagerMap.end() )
    {
        mDefaultDistributionType = defaultType;

        LAMA_LOG_INFO( logger, "set default communicator: " << mDefaultDistributionType )
    }
    else
    {
        setDefaultDistributionType( "NO" );

        LAMA_LOG_WARN( logger,
                       "default distribution " << defaultType << " not available, will be " << mDefaultDistributionType )
    }
}

/* ---- release  ------------------------------------------------------------------- */

void DistributionFactory::release()
{
    // release/free available context managers on the factory instance

    DistributionFactory& factory = getFactory();

    factory.mDistributionToManagerMap.clear();
}

/* ---- DistributionFactory constructor/destructor --------------------------------- */

DistributionFactory::DistributionFactory()
{
    mDefaultDistributionType = "";
}

DistributionFactory::~DistributionFactory()
{
    LAMA_LOG_INFO( logger, "~DistributionFactory" )
}

/* ---- DistributionFactory methods  ----------------------------------------------- */

void DistributionFactory::addDistributionManager(
    const std::string& type,
    boost::shared_ptr<DistributionManager> commManager )
{
    // be careful about logging, method might be called during static initialization

    DistributionToManagerMap::const_iterator it = mDistributionToManagerMap.find( type );

    if ( it != mDistributionToManagerMap.end() )
    {
        // DistributionManager might be replaced, e.g. for Host context (DefaultHost or CUDAHost)

        LAMA_LOG_INFO( logger, "Distribution manager replaced for communicator type " << type )
    }

    mDistributionToManagerMap[type] = commManager; // ownership of pointer is taken over
}

}
