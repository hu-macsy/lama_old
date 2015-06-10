/**
 * @file RegionEntry.hpp
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
 * @brief Definition of struct that contains all relevant properties of a region
 * @author Thomas Brandes
 * @date 11.06.2015
 */

#pragma once

#include <string>

namespace tracing
{

typedef unsigned int VTRegionId;

/** A region contains information about it source code location and about counters
 *  for time spent in it ( inclusive + exclusive)
 */

struct RegionEntry
{
    RegionEntry()
    {
        mLastTime = 0.0;
        mInclusiveTime = 0.0;
        mExclusiveTime = 0.0;
        mCalls = 0;

        mVTId = 0;
    }

    const char* getRegionName() const
    {
        return mName.c_str();
    }
    const char* getFileName() const
    {
        return mFile;
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

    /** subtract time of called regions to get exlusive time. */

    void subRegionCall( double spentTime )
    {
        mExclusiveTime -= spentTime;
    }

    std::string mName; // name of the region

    int mLine; // line where region starts

    VTRegionId mVTId; // region number given by VampirTrace

    const char* mFile; // file name where region is defined

private:

    int mCalls; //!< count number of calls for a region

    double mLastTime; //!< time of last call
    double mInclusiveTime; //!< time totally spent in a region
    double mExclusiveTime; //!< time exclusively spent in a region
};

}
