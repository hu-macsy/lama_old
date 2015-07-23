/**
 * @file CommunicatorFactory.cpp
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
 * @brief CommunicatorFactory.cpp
 * @author Thomas Brandes, Jiri Kraus
 * @date 23.02.2011
 * @since 1.0.0
 */

// hpp
#include <lama/CommunicatorFactory.hpp>

// others
#include <lama/CommunicatorManager.hpp>
#include <lama/exception/Exception.hpp>

// logging
#include <logging/logging.hpp>

// boost
#include <boost/shared_ptr.hpp>

using namespace std;
using namespace boost;

/** String for the default communicator. */

#define LAMA_ENV_FOR_COMMUNICATOR "LAMA_COMMUNICATOR"

namespace lama
{

class Communicator;

/* ---- Static variables of CommunicatorFactory  ----------------------------------- */

LAMA_LOG_DEF_LOGGER( CommunicatorFactory::logger, "CommunicatorFactory" )

CommunicatorFactory& CommunicatorFactory::getFactory()
{
    // construct on first use idiom

    static std::auto_ptr<CommunicatorFactory> theCommunicatorFactory;

    // Note: using the auto_ptr guarantees that the instance will be freed at end of the program

    if( !theCommunicatorFactory.get() )
    {
        // printf("create the communicator factory instance\n");

        theCommunicatorFactory = std::auto_ptr<CommunicatorFactory>( new CommunicatorFactory() );
    }

    // printf("get the factory\n");

    return *theCommunicatorFactory;
}

CommunicatorPtr CommunicatorFactory::get( const std::string& type, int& argc, char** & argv )
{
    LAMA_LOG_TRACE( logger, "Get communicator of type " << type )

    const CommunicatorFactory& factory = getFactory();

    CommunicatorToManagerMap::const_iterator id = factory.mCommunicatorToManagerMap.find( type );

    if( id != factory.mCommunicatorToManagerMap.end() )
    {
        return id->second->getCommunicator( argc, argv );
    }

    const std::string& defaultType = factory.getDefaultCommunicatorType();

    LAMA_LOG_WARN( logger, "Communicator for " << type << " not supported, will use " << defaultType )

    // try to get a default communicator

    id = factory.mCommunicatorToManagerMap.find( factory.getDefaultCommunicatorType() );

    if( id != factory.mCommunicatorToManagerMap.end() )
    {
        return id->second->getCommunicator( argc, argv );
    }

    COMMON_THROWEXCEPTION( "Communicator for " << type << " not supported, no manager registered" )
}

CommunicatorPtr CommunicatorFactory::get( const std::string& type )
{
    // provide NULL arguments

    int argc = 0;
    char** argv = NULL;
    return get( type, argc, argv );
}

CommunicatorPtr CommunicatorFactory::get()
{
    return get( getFactory().getDefaultCommunicatorType() );
}

const std::string& CommunicatorFactory::getDefaultCommunicatorType() const
{
    if( mDefaultCommunicatorType.length() == 0 )
    {
        // default communicator has not been defined yet, so do it

        if( getenv( LAMA_ENV_FOR_COMMUNICATOR ) )
        {
            setDefaultCommunicatorType( getenv( LAMA_ENV_FOR_COMMUNICATOR ) );
        }
        else
        {
            setDefaultCommunicatorType();
        }
    }

    return mDefaultCommunicatorType;
}

boost::shared_ptr<CommunicatorManager> CommunicatorFactory::getCommunicatorManager( const std::string& type )
{
    LAMA_ASSERT( mCommunicatorToManagerMap.find( type ) != mCommunicatorToManagerMap.end(),
                 "Communicator for " << type << " not supported, no manager registered" )
    return mCommunicatorToManagerMap[type];
}

void CommunicatorFactory::setDefaultCommunicatorType( const std::string& defaultType ) const
{
    if( mCommunicatorToManagerMap.find( defaultType ) != mCommunicatorToManagerMap.end() )
    {
        mDefaultCommunicatorType = defaultType;

        LAMA_LOG_INFO( logger, "set default communicator: " << mDefaultCommunicatorType )
    }
    else
    {
        setDefaultCommunicatorType();

        LAMA_LOG_WARN( logger,
                       "default communicator " << defaultType << " not available, will be " << mDefaultCommunicatorType )
    }
}

void CommunicatorFactory::setDefaultCommunicatorType() const
{
    if( mCommunicatorToManagerMap.find( "MPI" ) != mCommunicatorToManagerMap.end() )
    {
        mDefaultCommunicatorType = "MPI";
    }
    else if( mCommunicatorToManagerMap.find( "none" ) != mCommunicatorToManagerMap.end() )
    {
        mDefaultCommunicatorType = "none";
    }
    else
    {
        COMMON_THROWEXCEPTION( "No default communicator available" )
    }
}

/* ---- release  ------------------------------------------------------------------- */

void CommunicatorFactory::release()
{
    // release/free available context managers on the factory instance

    CommunicatorFactory& factory = getFactory();

    factory.mCommunicatorToManagerMap.clear();
}

/* ---- CommunicatorFactory constructor/destructor --------------------------------- */

CommunicatorFactory::CommunicatorFactory()
{
}

CommunicatorFactory::~CommunicatorFactory()
{
    LAMA_LOG_INFO( logger, "~CommunicatorFactory" )
}

/* ---- CommunicatorFactory methods  ----------------------------------------------- */

void CommunicatorFactory::addCommunicatorManager(
    const std::string& type,
    boost::shared_ptr<CommunicatorManager> commManager )
{
    // be careful about logging, method might be called during static initialization

    CommunicatorToManagerMap::const_iterator it = mCommunicatorToManagerMap.find( type );

    if( it != mCommunicatorToManagerMap.end() )
    {
        // CommunicatorManager might be replaced, e.g. for Host context (DefaultHost or CUDAHost)

        LAMA_LOG_INFO( logger, "Communicator manager replaced for communicator type " << type )
    }

    mCommunicatorToManagerMap[type] = commManager; // ownership of pointer is taken over
}

}
