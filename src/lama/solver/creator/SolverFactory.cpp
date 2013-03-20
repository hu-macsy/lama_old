/**
 * @file SolverFactory.cpp
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
 * @brief SolverFactory.cpp
 * @author Kai Buschulte
 * @date 06.06.2012
 * $Id$
 */

// hpp
#include <lama/solver/creator/SolverFactory.hpp>

#include <lama/exception/Exception.hpp>
#include <lama/exception/LAMAAssert.hpp>

namespace lama
{

LAMA_LOG_DEF_LOGGER( SolverFactory::logger, "Solver.SolverFactory" );

SolverFactory::SolverFactory()
{
}

SolverFactory::~SolverFactory()
{
    LAMA_LOG_DEBUG( logger, "~SolverFactory" );
}

void SolverFactory::addSolverCreator( const std::string& type, SolverCreator::RuleType& rule )
{
    if( mCreatorRuleMap.find( type ) )
    {
        LAMA_LOG_WARN( logger, "SolverCreator " << type << " wants to register twice." );
        return;
    }

    LAMA_LOG_DEBUG( logger, "Register rule of SolverCreator " << type << "." );

    mCreatorRuleMap.add( type.c_str(), rule );

    LAMA_LOG_DEBUG( logger, "SolverCreator type " << type << " registered." );
}

SolverFactory& SolverFactory::getFactory()
{
    static std::auto_ptr<SolverFactory> factory;

    if( !factory.get() )
    {
        factory = std::auto_ptr<SolverFactory>( new SolverFactory() );
    }
    return *factory;
}

const SolverFactory::TypeToCreatorMap& SolverFactory::getCreatorRuleSymbols()
{
    LAMA_LOG_DEBUG( logger, "Returning rule symbols map " << mCreatorRuleMap << "." );
    return mCreatorRuleMap;
}

void SolverFactory::addSolver( Solver* solver )
{
    SolverPtr solverPtr( solver );
    LAMA_ASSERT( solverPtr, "Solver not set." );

    const std::string& name = solverPtr->getId();

    if( hasSolver( name ) )
    {
        LAMA_LOG_WARN( logger, "Solver " << name << " wants to register with an existing name" );
        return;
    }

    mSolverInstanceMap.add( name, solverPtr );

    LAMA_LOG_INFO( logger, "Solver " << name << " registered " );
}

bool SolverFactory::hasSolver( const std::string& solverName ) const
{
    if( mSolverInstanceMap.find( solverName ) != NULL )
    {
        return true;
    }

    return false;
}

SolverPtr SolverFactory::getSolver( const std::string& solverName )
{
    if( !hasSolver( solverName ) )
    {
        LAMA_THROWEXCEPTION( "No solver instance named " << solverName << " registered." );
    }

    return mSolverInstanceMap.at( solverName );
}

const SolverFactory::SolverInstanceMap& SolverFactory::getSolverInstanceMap()
{
    return mSolverInstanceMap;
}

/* ---- release  ------------------------------------------------------------------- */

void SolverFactory::release()
{
    // release/free available context managers on the factory instance

    SolverFactory& factory = getFactory();

    LAMA_LOG_DEBUG( logger, "release creator rules and solver instances" );

    factory.mCreatorRuleMap.clear();
    factory.mSolverInstanceMap.clear();

    LAMA_LOG_DEBUG( logger, "released all creator rules and solver instances" );
}

}
