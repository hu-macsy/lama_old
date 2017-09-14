/**
 * @file BenchmarkRegistry.cpp
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
 * @brief BenchmarkRegistry.cpp
 * @author jiri
 * @date 06.04.2011
 */
/**
 * @file BenchmarkRegistry.cpp
 * @author jiri
 * Created on: 04.05.2010
 */

#include <scai/benchmark/BenchmarkRegistry.hpp>
#include <scai/benchmark/BFError.hpp>
#include <scai/benchmark/frame_stdlib.hpp>

#include <iostream>
#include <sstream>

namespace bf
{

BenchmarkRegistry* BenchmarkRegistry::mInstance = 0;

bool BenchmarkRegistry::isFreed = false;

AbstractBenchmarkCreator::~AbstractBenchmarkCreator()
{
}

BenchmarkRegistry::~BenchmarkRegistry()
{
}

BenchmarkRegistry& BenchmarkRegistry::getRegistry()
{
    static CGuard g;
    if( isFreed )
    {
        throw BFError( "BenchmarkRegistry is already freed. You can't get an instance." );
    }
    if( mInstance == 0 )
    {
        mInstance = new BenchmarkRegistry();
    }
    return *mInstance;
}

void BenchmarkRegistry::freeRegistry()
{
    if( mInstance != 0 )
    {
        delete mInstance;
        mInstance = 0;
        isFreed = true;
    }
}

AbstractBenchmarkCreator* BenchmarkRegistry::getCreatorById( const std::string& id ) const
{
    std::map<std::string,AbstractBenchmarkCreator*>::const_iterator loc = mCreatorMap.find( id );
    if( loc == mCreatorMap.end() )
    {
        throw BFError( "No creator Function for Benchmark " + id + " available" );
    }
    return loc->second;
}

std::auto_ptr<Benchmark> BenchmarkRegistry::createBenchmark( const std::string& id ) const
{
    std::string name = "";
    std::string args = "";

    std::string::size_type openBraked = id.find( '(', 0 );
    if( openBraked != std::string::npos )
    {
        ++openBraked;
        std::string::size_type closedBraked = id.find( ')', openBraked );
        if( closedBraked != std::string::npos )
        {
            name = id.substr( 0, openBraked - 1 );
            args = id.substr( openBraked, closedBraked - openBraked );
        }
        else
        {
            std::stringstream message;
            message << "BenchmarkRegistry::createBenchmark: missing ')' in " << "Benchmark name " << id;
            throw BFError( message.str() );
        }
        trimm( name );
        trimm( args );
        if( name.empty() )
        {
            std::stringstream message;
            message << "BenchmarkRegistry::createBenchmark: missing name in " << "Benchmark name " << id;
            throw BFError( message.str() );
        }
        if( args.empty() )
        {
            std::stringstream message;
            message << "BenchmarkRegistry::createBenchmark: missing arguments " << "in Benchmark name " << id;
            throw BFError( message.str() );
        }
    }
    else if( id.find( ')', 0 ) != std::string::npos )
    {
        std::stringstream message;
        message << "BenchmarkRegistry::createBenchmark: missing '(' in Benchmark " << "name " << id;
        throw BFError( message.str() );
    }
    else
    {
        name = id;
    }
    AbstractBenchmarkCreator* bc = getCreatorById( name );

    if( args.empty() )
    {
        return std::auto_ptr<Benchmark>( bc->create() );
    }
    else
    {
        return std::auto_ptr<Benchmark>( bc->create( args ) );
    }
}

std::auto_ptr<Benchmark> BenchmarkRegistry::createBenchmark( const std::string& id, const std::string& args ) const
{
    AbstractBenchmarkCreator* bc = getCreatorById( id );
    return std::auto_ptr<Benchmark>( bc->create( args ) );
}

void BenchmarkRegistry::destroyBenchmark( Benchmark* bench ) const
{
    delete bench;
}

void BenchmarkRegistry::add( const std::string& id, AbstractBenchmarkCreator* creator )
{
    mCreatorMap.insert( std::make_pair( id, creator ) );
}

bool BenchmarkRegistry::has( const std::string& id ) const
{
    std::string name = "";

    std::string::size_type openBraked = id.find( '(', 0 );
    if( openBraked != std::string::npos )
    {
        std::string::size_type closedBraked = id.find( ')', openBraked );
        if( closedBraked != std::string::npos )
        {
            name = id.substr( 0, openBraked );
        }
        else
        {
            std::stringstream message;
            message << "BenchmarkRegistry::createBenchmark: missing ')' in " << "Benchmark name " << id;
            throw BFError( message.str() );
        }
        trimm( name );
        if( name.empty() )
        {
            std::stringstream message;
            message << "BenchmarkRegistry::createBenchmark: missing name in " << "Benchmark name " << id;
            throw BFError( message.str() );
        }
    }
    else if( id.find( ')', 0 ) != std::string::npos )
    {
        std::stringstream message;
        message << "BenchmarkRegistry::createBenchmark: missing '(' in Benchmark " << "name " << id;
        throw BFError( message.str() );
    }
    else
    {
        name = id;
    }
    std::map<std::string,AbstractBenchmarkCreator*>::const_iterator loc = mCreatorMap.find( name );
    return loc != mCreatorMap.end();
}

BenchmarkRegistry::BenchmarkRegistry()
{

}

BenchmarkRegistry::CGuard::CGuard()
{
}

BenchmarkRegistry::CGuard::~CGuard()
{
    if( BenchmarkRegistry::mInstance != 0 )
    {
        delete BenchmarkRegistry::mInstance;
        BenchmarkRegistry::mInstance = 0;
    }
}

BenchmarkRegistry::const_iterator BenchmarkRegistry::begin() const
{
    return mCreatorMap.begin();
}

BenchmarkRegistry::const_iterator BenchmarkRegistry::end() const
{
    return mCreatorMap.end();
}
} // namespace bf
