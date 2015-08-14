/**
 * @file AMGSetupFactory.cpp
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
 * @brief AMGSetupFactory.cpp
 * @author Thomas Brandes, Jiri Kraus
 * @date 23.02.2011
 * @since 1.0.0
 */

// hpp
#include <scai/lama/solver/AMGSetupFactory.hpp>

// others
#include <scai/lama/solver/AMGSetupManager.hpp>

#include <scai/lama/exception/Exception.hpp>

// logging
#include <scai/logging.hpp>

// boost
#include <scai/common/shared_ptr.hpp>

using namespace std;
using namespace common;

/** String for the default AMG Setup. */

#define LAMA_ENV_FOR_COMMUNICATOR "LAMA_COMMUNICATOR"

namespace scai
{

namespace lama
{

class AMGSetup;

/* ---- Static variables of AMGSetupFactory  ----------------------------------- */

SCAI_LOG_DEF_LOGGER( AMGSetupFactory::logger, "AMGSetupFactory" )

AMGSetupFactory& AMGSetupFactory::getFactory()
{
    // construct on first use idiom

    static common::unique_ptr<AMGSetupFactory> theAMGSetupFactory;

    // Note: using the auto_ptr guarantees that the instance will be freed at end of the program

    if( !theAMGSetupFactory.get() )
    {
        // printf("create the AMGSetup factory instance\n");

        theAMGSetupFactory.reset( new AMGSetupFactory() );
    }

    // printf("get the factory\n");

    return *theAMGSetupFactory;
}

AMGSetupPtr AMGSetupFactory::get( const std::string& type )
{
    SCAI_LOG_TRACE( logger, "Get AMGSetup of type " << type )

    const AMGSetupFactory& factory = getFactory();

    AMGSetupToManagerMap::const_iterator id = factory.mAMGSetupToManagerMap.find( type );

    if( id != factory.mAMGSetupToManagerMap.end() )
    {
        return id->second->getAMGSetup();
    }

    const std::string& defaultType = factory.getDefaultAMGSetupType();

    SCAI_LOG_WARN( logger, "AMGSetup for " << type << " not supported, will use " << defaultType )

    // try to get a default AMGSetup

    id = factory.mAMGSetupToManagerMap.find( factory.getDefaultAMGSetupType() );

    if( id != factory.mAMGSetupToManagerMap.end() )
    {
        return id->second->getAMGSetup();
    }

    COMMON_THROWEXCEPTION( "AMGSetup for " << type << " not supported, no manager registered" )
}

AMGSetupPtr AMGSetupFactory::get()
{
    return get( getFactory().getDefaultAMGSetupType() );
}

const std::string& AMGSetupFactory::getDefaultAMGSetupType() const
{
    if( mDefaultAMGSetupType.length() == 0 )
    {
        // default AMGSetup has not been defined yet, so do it

        if( getenv( LAMA_ENV_FOR_COMMUNICATOR ) )
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

common::shared_ptr<AMGSetupManager> AMGSetupFactory::getAMGSetupManager( const std::string& type )
{
    SCAI_ASSERT( mAMGSetupToManagerMap.find( type ) != mAMGSetupToManagerMap.end(),
                 "AMGSetup for " << type << " not supported, no manager registered" )
    return mAMGSetupToManagerMap[type];
}

void AMGSetupFactory::setDefaultAMGSetupType( const std::string& defaultType ) const
{
    if( mAMGSetupToManagerMap.find( defaultType ) != mAMGSetupToManagerMap.end() )
    {
        mDefaultAMGSetupType = defaultType;

        SCAI_LOG_INFO( logger, "set default AMGSetup: " << mDefaultAMGSetupType )
    }
    else
    {
        setDefaultAMGSetupType();

        SCAI_LOG_WARN( logger,
                       "default AMGSetup " << defaultType << " not available, will be " << mDefaultAMGSetupType )
    }
}

void AMGSetupFactory::setDefaultAMGSetupType() const
{
    if( mAMGSetupToManagerMap.find( "MPI" ) != mAMGSetupToManagerMap.end() )
    {
        mDefaultAMGSetupType = "MPI";
    }
    else if( mAMGSetupToManagerMap.find( "none" ) != mAMGSetupToManagerMap.end() )
    {
        mDefaultAMGSetupType = "none";
    }
    else
    {
        COMMON_THROWEXCEPTION( "No default AMGSetup available" )
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
    SCAI_LOG_INFO( logger, "~AMGSetupFactory" )
}

/* ---- AMGSetupFactory methods  ----------------------------------------------- */

void AMGSetupFactory::addAMGSetupManager( const std::string& type, common::shared_ptr<AMGSetupManager> commManager )
{
    // be careful about logging, method might be called during static initialization

    AMGSetupToManagerMap::const_iterator it = mAMGSetupToManagerMap.find( type );

    if( it != mAMGSetupToManagerMap.end() )
    {
        // AMGSetupManager might be replaced, e.g. for Host context (DefaultHost or CUDAHost)

        SCAI_LOG_INFO( logger, "AMGSetup manager replaced for AMGSetup type " << type )
    }

    mAMGSetupToManagerMap[type] = commManager; // ownership of pointer is taken over
}

} /* end namespace lama */

} /* end namespace scai */
