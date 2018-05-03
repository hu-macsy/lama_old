/**
 * @file LamaTiming.hpp
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
 * @brief Class that is useful for timing in MPI programs
 * @author Thomas Brandes
 * @date 09.07.2013
 */

#pragma once

#include <scai/common/Walltime.hpp>
#include <iostream>
#include <scai/dmemo/Communicator.hpp>

/** Class for Timing.
 *
 *  \code
 *      LamaTiming time( comm, "Redistribution" )
 *  \endcode
 */

class LamaTiming
{

public:

    /** Constructor. */

    LamaTiming( const scai::dmemo::Communicator& comm, const char* name );

    /** Destructor, prints timing on root processor */

    ~LamaTiming();

    /** get wall time spent after calling the constructor */

    double getTime() const;

private:

    const scai::dmemo::Communicator& mComm;
    const char* mName;
    double mStart;
};

/* ---------------------------------------------------------------------------- */

LamaTiming::LamaTiming( const scai::dmemo::Communicator& comm, const char* name ) :
    mComm( comm ),
    mName( name )
{
    mStart = scai::common::Walltime::get();
}

double LamaTiming::getTime() const
{
    double myTime = scai::common::Walltime::get() - mStart;
    // can be that max is not available if double is not supported
    double maxTime = static_cast<double>( mComm.max( scai::DefaultReal( myTime ) ) );

    return maxTime;
}

LamaTiming::~LamaTiming()
{
    double time = getTime();

    if ( mComm.getRank() == 0 )
    {
        std::cout << mName << ": took " << time << " seconds" << std::endl;
    }
}
