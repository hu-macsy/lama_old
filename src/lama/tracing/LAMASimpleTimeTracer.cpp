/**
 * @file LAMASimpleTimeTracer.cpp
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
 * @brief LAMASimpleTimeTracer.cpp
 * @author Lauretta Schubert
 * @date 08.11.2011
 * $Id$
 */

// hpp
#include <lama/tracing/LAMASimpleTimeTracer.hpp>

#include <iostream>

std::list<std::pair<std::string,double> > LAMASimpleTimeTracer::timerList;
boost::mutex LAMASimpleTimeTracer::access_mutex;

double LAMASimpleTimeTracer::spentLast( const char* name )
{
    boost::mutex::scoped_lock scoped_lock( access_mutex );

    const std::string aggregated = "_aggregated";
    const std::string nameToSearch = name;
    const std::string aggregatedName = name + aggregated;

    const int suffixLength = aggregated.length();
    std::list<std::pair<std::string,double> >::iterator foundIndex;
    bool found = false;

    //search aggregated name
    std::list<std::pair<std::string,double> >::iterator end = timerList.end();
    for ( std::list<std::pair<std::string,double> >::iterator it = timerList.begin(); it != end; ++it )
    {
        if ( it->first.length() > (unsigned int) suffixLength
                && it->first.substr( it->first.length() - suffixLength ) != aggregated )
        {
            break;
        }
        if ( it->first == aggregatedName )
        {
            found = true;
            foundIndex = it;
            break;
        }
    }

    if ( !found )
    {
        std::pair<std::string,double> newEntry = std::pair<std::string,double>( aggregatedName, 0.0 );
        timerList.push_front( newEntry );
        foundIndex = timerList.begin();
    }

    end = timerList.end();
    double spentLastTime = 0.0;
    for ( std::list<std::pair<std::string,double> >::iterator it = timerList.begin(); it != end;
            /*it is incremented in loop body*/)
    {
        if ( it->first == nameToSearch )
        {
            spentLastTime += it->second;

            // first increment iteration, remove afterwards the element in list
            std::list<std::pair<std::string,double> >::iterator elementToRemove = it;
            ++it;
            timerList.erase( elementToRemove );
        }
        else
        {
            ++it;
        }
    }
    foundIndex->second += spentLastTime;

    return spentLastTime;
}

void LAMASimpleTimeTracer::printTimer()
{
    //TODO implement
}
