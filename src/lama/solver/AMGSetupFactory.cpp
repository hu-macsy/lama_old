/**
 * @file AMGSetupFactory.cpp
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
 * @brief AMGSetupFactory.cpp
 * @author Thomas Brandes, Jiri Kraus
 * @date 23.02.2011
 * $Id$
 */

// hpp
#include <lama/solver/AMGSetupFactory.hpp>

// others
#include <lama/solver/AMGSetupManager.hpp>

#include <lama/exception/Exception.hpp>

// logging
#include <logging/logging.hpp>

// boost
#include <boost/shared_ptr.hpp>

using namespace std;
using namespace boost;

/** String for the default AMG Setup. */

#define LAMA_ENV_FOR_COMMUNICATOR "LAMA_COMMUNICATOR"

namespace lama
{

class AMGSetup;

/* ---- Static variables of AMGSetupFactory  ----------------------------------- */

LAMA_LOG_DEF_LOGGER( AMGSetupFactory::logger, "AMGSetupFactory" )

AMGSetupFactory& AMGSetupFactory::getFactory()
{
    // construct on first use idiom

    static std::auto_ptr<AMGSetupFactory> theAMGSetupFactory;

    // Note: using the auto_ptr guarantees that the instance will be freed at end of the program

    if ( !theAMGSetupFactory.get() )
    {
        // printf("create the AMGSetup factory instance\n");

        theAMGSetupFactory = std::auto_ptr<AMGSetupFactory>( new AMGSetupFactory() );
    }

    // printf("get the factory\n");

    return *theAMGSetupFactory;
}

AMGSetupPtr AMGSetupFactory::get( const std::string& type )
{
    LAMA_LOG_TRACE( logger, "Get AMGSetup of type " << type )

    const AMGSetupFactory& factory = getFactory();

    AMGSetupToManagerMap::const_iterator id = factory.mAMGSetupToManagerMap.find( type );

    if ( id != factory.mAMGSetupToManagerMap.end() )
    {
        return id->second->getAMGSetup();
    }

    const std::string& defaultType = factory.getDefaultAMGSetupType();

    LAMA_LOG_WARN( logger, "AMGSetup for " << type << " not supported, will use " << defaultType )

    // try to get a default AMGSetup

    id = factory.mAMGSetupToManagerMap.find( factory.getDefaultAMGSetupType() );

    if ( id != factory.mAMGSetupToManagerMap.end() )
    {
        return id->second->getAMGSetup();
    }

    LAMA_THROWEXCEPTION( "AMGSetup for " << type << " not supported, no manager registered" )
}

AMGSetupPtr AMGSetupFactory::get()
{
    return get( getFactory().getDefaultAMGSetupType() );
}

const std::string& AMGSetupFactory::getDefaultAMGSetupType() const
{
    if ( mDefaultAMGSetupType.length() == 0 )
    {
        // default AMGSetup has not been defined yet, so do it

        if ( getenv( LAMA_ENV_FOR_COMMUNICATOR ) )
        {
            setDefaultAMGSetupType( getenv( LAMA_ENV_FOR_COMMUNICATOR ) );
        }
        else
        {
            setDefaultAMGSetupType();
        }
    }

    return mDefaultAMGSetupType;
}

boost::shared_ptr<AMGSetupManager> AMGSetupFactory::getAMGSetupManager( const std::string& type )
{
    LAMA_ASSERT( mAMGSetupToManagerMap.find( type ) != mAMGSetupToManagerMap.end(),
                 "AMGSetup for " << type << " not supported, no manager registered" )
    return mAMGSetupToManagerMap[type];
}

void AMGSetupFactory::setDefaultAMGSetupType( const std::string& defaultType ) const
{
    if ( mAMGSetupToManagerMap.find( defaultType ) != mAMGSetupToManagerMap.end() )
    {
        mDefaultAMGSetupType = defaultType;

        LAMA_LOG_INFO( logger, "set default AMGSetup: " << mDefaultAMGSetupType )
    }
    else
    {
        setDefaultAMGSetupType();

        LAMA_LOG_WARN( logger,
                       "default AMGSetup " << defaultType << " not available, will be " << mDefaultAMGSetupType )
    }
}

void AMGSetupFactory::setDefaultAMGSetupType() const
{
    if ( mAMGSetupToManagerMap.find( "MPI" ) != mAMGSetupToManagerMap.end() )
    {
        mDefaultAMGSetupType = "MPI";
    }
    else if ( mAMGSetupToManagerMap.find( "none" ) != mAMGSetupToManagerMap.end() )
    {
        mDefaultAMGSetupType = "none";
    }
    else
    {
        LAMA_THROWEXCEPTION( "No default AMGSetup available" )
    }
}

/* ---- release  ------------------------------------------------------------------- */

void AMGSetupFactory::release()
{
    // release/free available context managers on the factory instance

    AMGSetupFactory& factory = getFactory();

    factory.mAMGSetupToManagerMap.clear();
}

/* ---- AMGSetupFactory constructor/destructor --------------------------------- */

AMGSetupFactory::AMGSetupFactory()
{
}

AMGSetupFactory::~AMGSetupFactory()
{
    LAMA_LOG_INFO( logger, "~AMGSetupFactory" )
}

/* ---- AMGSetupFactory methods  ----------------------------------------------- */

void AMGSetupFactory::addAMGSetupManager( const std::string& type, boost::shared_ptr<AMGSetupManager> commManager )
{
    // be careful about logging, method might be called during static initialization

    AMGSetupToManagerMap::const_iterator it = mAMGSetupToManagerMap.find( type );

    if ( it != mAMGSetupToManagerMap.end() )
    {
        // AMGSetupManager might be replaced, e.g. for Host context (DefaultHost or CUDAHost)

        LAMA_LOG_INFO( logger, "AMGSetup manager replaced for AMGSetup type " << type )
    }

    mAMGSetupToManagerMap[type] = commManager; // ownership of pointer is taken over
}

}
