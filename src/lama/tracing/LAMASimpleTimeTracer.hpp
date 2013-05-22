/**
 * @file LAMASimpleTimeTracer.hpp
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
 * @brief LAMASimpleTimeTracer.hpp
 * @author Lauretta Schubert
 * @date 08.11.2011
 * $Id$
 */
#ifndef LAMA_LAMASIMPLETIMETRACER_HPP_
#define LAMA_LAMASIMPLETIMETRACER_HPP_

// for dll_import
#include <lama/config.hpp>

// base classes
#include <lama/tracing/LAMABaseTracer.hpp>

// boost
#include <boost/thread/mutex.hpp>

#include <list>
#include <string>
#include <omp.h>

/**
 * @brief lama::LAMASimpleTimeTracer
 */

class LAMASimpleTimeTracer: public LAMABaseTracer
{
public:

    inline LAMASimpleTimeTracer( const char* name, const char* file, int lno );

    inline ~LAMASimpleTimeTracer();

    /* be careful this is the spentLastTime since the last call of this function, */
    /* not only the one of the last function call */
    static double spentLast( const char* name );

    static void printTimer();

private:

    static std::list<std::pair<std::string,double> > timerList;

    static boost::mutex access_mutex; // needed to make time thread-safe

    const std::string mName;

    double mStopTime;

    const double mStartTime; // mStartTime should be initialized last
};

inline LAMASimpleTimeTracer::LAMASimpleTimeTracer( const char* name, const char* /*file*/, int /*lno*/)
    : mName( name ), mStopTime( 0.0 ), mStartTime( omp_get_wtime() )
{
}

inline LAMASimpleTimeTracer::~LAMASimpleTimeTracer()
{
    mStopTime = omp_get_wtime();

    double runTime = mStopTime - mStartTime;
    if ( getRuntime() > 0.0 )
    {
        runTime = getRuntime();
    }

    std::pair<std::string,double> value = std::make_pair( mName, runTime );
    {
        boost::mutex::scoped_lock scoped_lock( access_mutex );
        timerList.push_back( value );
    }
}

#endif // LAMA_LAMASIMPLETIMETRACER_HPP_
