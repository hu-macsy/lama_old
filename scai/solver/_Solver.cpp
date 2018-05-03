/**
 * @file _Solver.cpp
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
 * @brief Implementation of methods for the base class _Solver.
 * @author Jiri Kraus
 * @date 08.06.2011
 */

// hpp
#include <scai/solver/_Solver.hpp>

// local library
#include <scai/solver/logger/CommonLogger.hpp>
#include <scai/solver/logger/Timer.hpp>

namespace scai
{

namespace solver
{

/* ========================================================================= */
/*    Implementation of methods for _Solver                                  */
/* ========================================================================= */

SCAI_LOG_DEF_LOGGER( _Solver::logger, "Solver" )

_Solver* _Solver::getSolver( const common::ScalarType scalarType, const std::string& solverType )
{
    return create( SolverCreateKeyType( scalarType, solverType ) );
}

_Solver::_Solver( const std::string& id ) : 

    mId( id ),

    mLogger( new CommonLogger( "dummyLog", LogLevel::noLogging,
                           LoggerWriteBehaviour::toConsoleOnly,
                           std::shared_ptr<Timer>( new Timer() ) ) )

{
    SCAI_LOG_INFO( _Solver::logger, "Solver id = " << mId << " created, dummy log" )
}

_Solver::_Solver( const std::string& id, LoggerPtr logger ) :

    mId( id ),
    mLogger( logger )
{
}

_Solver::_Solver( const _Solver& other ) : 

    mId( other.mId ),
    mLogger( other.mLogger )    // Be careful: no deep copy of the logger here
{
}

void _Solver::setId( const std::string& id ) 
{
    mId = id;
}

void _Solver::setLogger( LoggerPtr logger )
{
    mLogger = logger;
}

void _Solver::setLogLevel( LogLevel level )
{
    mLogger->setLogLevel( level );
}

void _Solver::getTypedCreateValues( std::vector<std::string>& values, const common::ScalarType stype )
{
    std::vector<SolverCreateKeyType> createValues;

    common::Factory<SolverCreateKeyType, _Solver*>::getCreateValues( createValues );  // all solvers ( valueType, solvertype )

    values.clear();

    for ( size_t i = 0; i < createValues.size(); ++i )
    {
        if ( createValues[i].first == stype )
        {
            // Solver for this value type
            values.push_back( createValues[i].second );
        }
    }
}

} /* end namespace solver */

} /* end namespace scai */
