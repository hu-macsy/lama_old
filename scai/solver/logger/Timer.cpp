/**
 * @file Timer.cpp
 *
 * @license
 * Copyright (c) 2009-2018
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the SCAI framework LAMA.
 *
 * LAMA is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * LAMA is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
 * more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with LAMA. If not, see <http://www.gnu.org/licenses/>.
 * @endlicense
 *
 * @brief Timer.cpp
 * @author Matthias Makulla
 * @date 06.04.2011
 */

// hpp
#include <scai/solver/logger/Timer.hpp>

// local library
#include <scai/common/macros/throw.hpp>
#include <scai/common/Walltime.hpp>

// std
#include <iostream>
#include <cstdio>
#include <sstream>

namespace scai
{

namespace solver
{

Timer::Timer()
{
}

Timer::~Timer()
{
}

void Timer::initialize( const std::string& timerId )
{
    MapIteratorType it = m_timerData.find( timerId );

    if ( it == m_timerData.end() )
    {
        m_timerData.insert( PairType( timerId, TimerData() ) );
    }
    else
    {
        COMMON_THROWEXCEPTION( "Tried to initialize an already existing timer. Timer ID: " << timerId );
    }
}

void Timer::start( const std::string& timerId )
{
    MapIteratorType it = m_timerData.find( timerId );

    if ( it == m_timerData.end() )
    {
        it = m_timerData.insert( PairType( timerId, TimerData() ) ).first;
    }

    TimerData& timer = it->second;

    if ( timer.isRunning )
    {
        COMMON_THROWEXCEPTION( "Tried to start an already started timer. Timer ID: " << timerId );
    }
    else
    {
        timer.startTime = common::Walltime::get();
        timer.isRunning = true;
    }
}

void Timer::stop( const std::string& timerId )
{
    MapIteratorType it = m_timerData.find( timerId );

    if ( it == m_timerData.end() )
    {
        COMMON_THROWEXCEPTION( "Tried to stop a nonexisting Timer. Timer ID: " << timerId );
    }

    TimerData& timer = it->second;

    if ( !( timer.isRunning ) )
    {
        COMMON_THROWEXCEPTION( "Tried to stop a not running Timer. Timer ID: " << timerId );
    }

    timer.isRunning = false;
    timer.totalTime = timer.totalTime + ( common::Walltime::get() - timer.startTime );
}

void Timer::reset( const std::string& timerId )
{
    MapIteratorType it = m_timerData.find( timerId );

    if ( it == m_timerData.end() )
    {
        COMMON_THROWEXCEPTION( "Tried to reset a nonexisting Timer. Timer ID: " << timerId );
    }

    TimerData& timer = it->second;
    timer.totalTime = 0.0;
    timer.startTime = common::Walltime::get();
}

double Timer::getTime( const std::string& timerId )
{
    MapIteratorType it = m_timerData.find( timerId );

    if ( it == m_timerData.end() )
    {
        COMMON_THROWEXCEPTION( "Tried to get time from a nonexisting Timer. Timer ID: " << timerId );
    }

    TimerData& timer = it->second;

    if ( !( timer.isRunning ) )
    {
        return timer.totalTime;
    }

    return timer.totalTime + ( common::Walltime::get() - timer.startTime );
}

void Timer::stopAndReset( const std::string& timerId )
{
    MapIteratorType it = m_timerData.find( timerId );

    if ( it == m_timerData.end() )
    {
        COMMON_THROWEXCEPTION( "Tried to stop and reset a nonexisting Timer. Timer ID: " << timerId );
    }

    TimerData& timer = it->second;
    timer.totalTime = 0.0;
    timer.startTime = 0.0;
    timer.isRunning = false;
}

} /* end namespace solver */

} /* end namespace scai */
