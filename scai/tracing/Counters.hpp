/**
 * @file Counters.hpp
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
 * @brief Definition of class for Performance Counters
 * @author Thomas Brandes
 * @date 23.06.2015
 */

#pragma once

// internal scai libraries
#include <scai/common/Walltime.hpp>
#include <scai/common/Printable.hpp>

// std
#include <ostream>

namespace scai
{

namespace tracing
{

static const int MAX_COUNTERS = 1;

static const int usedCounters = 1;    //!< number of counters used, <= MAX_COUNTERS

static inline int enabledCounters()
{
    return usedCounters;
}

class CounterArray : public scai::common::Printable
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

        return *this;
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
        return double( values[0] ) / double( common::Walltime::timerate() );
    }

    /** Return walltime as difference to previous counter values */

    double getWalltime( const CounterArray& other ) const
    {
        return double ( values[0] - other.values[0] ) / double( common::Walltime::timerate() );
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

} /* end namespace tracing */

} /* end namespace scai */
