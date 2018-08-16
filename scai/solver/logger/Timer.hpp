/**
 * @file Timer.hpp
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
 * @brief Class that provides a set of timers accessed by their ids.
 * @author Matthias Makulla
 * @date 06.04.2011
 */

#pragma once

// for dll_import
#include <scai/common/config.hpp>

// std
#include <map>
#include <string>
#include <memory>

namespace scai
{

namespace solver
{

class Timer;
typedef std::shared_ptr<Timer> TimerPtr;

/**
 * @brief Timer class that offers a set of timers identified by their ids.
 *
 * An object of this class provides routines to start and stop timers that
 * are identified by a string.
 */
class COMMON_DLL_IMPORTEXPORT Timer
{
public:

    /** Constructor of a new object for set of timers. */

    Timer();

    virtual ~Timer();

    /**
     * @brief Starts or resumes the timer.
     *
     * @param[in] timerId   the ID of the timer
     */
    void start( const std::string& timerId );

    /**
     * @brief Stops the timer and stores time measured time internally.
     *        Measurement may be resumed by calling start()
     *
     * @param[in] timerId   the ID of the timer
     * @throw Exception if timer has already been started
     */
    void stop( const std::string& timerId );

    /**
     * @brief Gets the elapsed time since the last reset
     *
     * @param[in] timerId   the ID of the timer
     * @return Time measured in seconds
     */
    double getTime( const std::string& timerId );

    /**
     * @brief Stops and resets the timer.
     *
     * @param[in] timerId   the ID of the timer
     *
     * Note: this routine might be typically called after a call of getTime
     */
    void stopAndReset( const std::string& timerId );

    /**
     * @brief Reset a timer.
     *
     * @param[in] timerId  the ID of an existing timer
     *
     * The timer might also be an already started timer.
     */
    void reset( const std::string& timerId );

private:

    void initialize( const std::string& timerId );

    struct TimerData
    {
        double startTime;
        double totalTime;
        bool isRunning;

        TimerData()
            : startTime( 0.0 ), totalTime( 0.0 ), isRunning( false )
        {
        }
    };

    typedef std::map<std::string, TimerData> MapType;
    typedef std::pair<std::string, TimerData> PairType;
    typedef MapType::iterator MapIteratorType;

    MapType m_timerData;

};

} /* end namespace solver */

} /* end namespace scai */
