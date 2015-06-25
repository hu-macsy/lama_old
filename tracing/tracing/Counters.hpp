/**
 * @file Counters.hpp
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
 * @brief Definition of class for Performance Counters
 * @author Thomas Brandes
 * @date 23.06.2015
 */

#pragma once

#include <common/Walltime.hpp>
#include <common/Printable.hpp>

#include <ostream>

namespace tracing
{
static const int MAX_COUNTERS = 1;

static const int usedCounters = 1;    //!< number of counters used, <= MAX_COUNTERS

static inline int enabledCounters()
{
    return usedCounters;
}

class CounterArray : public Printable
{

public:

    /** Constructor of counter sets the current values. */

    CounterArray( bool stampIt = false )
    {
        if ( stampIt )
        {
            stamp();
        }
        else
        {
            for ( int i = 0; i < usedCounters; ++i )
            {
                values[i] = 0;
            }
        }
    }

    ~CounterArray()
    {
    }

    CounterArray operator- ( const CounterArray& other ) const
    {
        CounterArray result;

        for ( int i = 0; i < usedCounters; ++i )
        {
            result.values[i] = this->values[i] - other.values[i];
        }

        return result;
    }

    /** Override default assignment operator, does not copy all */

    CounterArray& operator= ( const CounterArray& other )
    {
        for ( int i = 0; i < usedCounters; ++i )
        {
            values[i] = other.values[i];
        }
    }

    /** Override default copy constructor. */

    CounterArray( const CounterArray& other )
    {
        *this = other;   // uses assignment operator.
    }

    CounterArray& operator += ( const CounterArray& other )
    {
        for ( int i = 0; i < usedCounters; ++i )
        {
            values[i] += other.values[i];
        }

        return *this;
    }

    uint64_t operator[] ( int i ) const
    {
        return values[i];
    }

    double getWalltime() const
    {
        return values[0] / common::Walltime::timerate();
    }

    /** Return walltime as difference to previous counter values */

    double getWalltime( const CounterArray& other ) const
    {
        return ( values[0] - other.values[0] ) / common::Walltime::timerate();
    }

    /** This method writes just the counter values in a stream, separated by a given string. */

    void write( std::ostream& stream, const char* separator ) const
    {
        for ( int i = 0; i < usedCounters; ++i )
        {
            stream << values[i];

            if ( i + 1 < usedCounters )
            {
                stream << separator;
            }
        }
    }

    virtual void writeAt( std::ostream& stream ) const
    {
        stream << "CounterArray={ " ;
        write( stream, ", " );
        stream << " }";
    }

private:

    void stamp()
    {
        // currently we count only time ticks
        values[0] = common::Walltime::timestamp();
    }

    uint64_t values[MAX_COUNTERS];
};

}  // namespace
