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
 * @brief Timer implementation using OpenMP timing
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
 * @brief Timer implementation using OpenMP timing
 *        For further documentation see Timer interface
 */
class LAMA_DLL_IMPORTEXPORT OpenMPTimer: public Timer
{
public:

    OpenMPTimer();
    virtual ~OpenMPTimer();

    void initialize( const std::string& timerId );

    void start( const std::string& timerId );
    void stop( const std::string& timerId );

    void reset( const std::string& timerId );

    double getTime( const std::string& timerId );
    void stopAndReset( const std::string& timerId );

private:

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
