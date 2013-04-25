/**
 * @file OpenMPTimer.hpp
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
 * @brief Class that provides a set of timers accessed by their ids.
 * @author Matthias Makulla
 * @date 06.04.2011
 * $Id$
 */
#ifndef LAMA_OPENMPTIMER_HPP_
#define LAMA_OPENMPTIMER_HPP_

// for dll_import
#include <lama/config.hpp>

// base classes
#include <lama/solver/logger/Timer.hpp>

#include <map>

namespace lama
{

/**
 * @brief Timer class that offers a set of timers identified by their ids.
 *
 * An object of this class provides routines to start and stop timers that
 * are identified by a string.
 */
class LAMA_DLL_IMPORTEXPORT OpenMPTimer: public Timer
{
public:

    /** Constructor of a new object for set of timers. */

    OpenMPTimer();

    virtual ~OpenMPTimer();

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

    typedef std::map<std::string,TimerData> MapType;
    typedef std::pair<std::string,TimerData> PairType;
    typedef MapType::iterator MapIteratorType;

    MapType m_timerData;

};

} // namespace lama

#endif // LAMA_OPENMPTIMER_HPP_
