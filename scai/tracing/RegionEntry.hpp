/**
 * @file RegionEntry.hpp
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
 * @brief Definition of struct that contains all relevant properties of a region
 * @author Thomas Brandes
 * @date 11.06.2015
 */

#pragma once

// std
#include <string>
#include <fstream>
#include <iostream>
#include <iomanip>

namespace scai
{

namespace tracing
{

typedef unsigned int VTRegionId;

/** A region contains information about it source code location and about counters
 *  for time spent in it ( inclusive + exclusive)
 */

class RegionEntry
{
public:

    RegionEntry()
    {
        mLastTime = 0.0;
        mInclusiveTime = 0.0;
        mExclusiveTime = 0.0;
        mCalls = 0;
        mVTId = 0;
        mFirst = true;    // for first access
    }

    ~RegionEntry()
    {
    }

    const char* getRegionName() const
    {
        return mName.c_str();
    }

    const char* getFileName() const
    {
        return mFile;
    }

    int getFileToken() const
    {
        return 0;
    }

    int getLine() const
    {
        return mLine;
    }

    int getCalls() const
    {
        return mCalls;
    }

    double getLastTime() const
    {
        return mLastTime;
    }

    double getInclusiveTime() const
    {
        return mInclusiveTime;
    }

    double getExclusiveTime() const
    {
        return mExclusiveTime;
    }

    void addCall( double spentTime )
    {
        mLastTime = spentTime;
        mInclusiveTime += spentTime;
        mExclusiveTime += spentTime;
        mCalls++;
    }

    bool firstAccess()
    {
        bool is = mFirst;
        mFirst = false;
        return is;
    }

    /** subtract time of called regions to get exlusive time. */

    void subRegionCall( double spentTime )
    {
        mExclusiveTime -= spentTime;
    }

    std::string mName; // name of the region

    int mLine; // line where region starts

    VTRegionId mVTId; // region number given by VampirTrace

    const char* mFile; // file name where region is defined

    void printTime( std::ostream& outfile ) const
    {
        outfile << "Time " << mName << " (in ms) : ";
        outfile << "#calls = " << mCalls;
        outfile << ", inclusive = " << std::fixed << std::setprecision( 6 ) << ( mInclusiveTime * 1000.0 );
        outfile << ", exclusive = " << std::fixed << std::setprecision( 6 ) << ( mExclusiveTime * 1000.0 );
        outfile << std::endl;
    }

    /*
    virtual void writeAt( std::ostream& outfile ) const
    {
        outfile << "Region( name = " << mName << " )";
    }
    */

private:

    int mCalls;            //!< counts number of calls for this region

    double mLastTime;      //!< time spent on latest call

    double mInclusiveTime; //!< time totally spent in a region
    double mExclusiveTime; //!< time exclusively spent in a region

    bool mFirst;           //!< only true for first access
};

} /* end namespace tracing */

} /* end namespace scai */
