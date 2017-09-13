/**
 * @file Benchmark.cpp
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
 * @brief Benchmark.cpp
 * @author jiri
 * @date 06.04.2011
 * $Id$
 */
/**
 * @file Benchmark.cpp
 * @author jiri
 * Created on: 04.05.2010
 */
#include <scai/benchmark/Benchmark.hpp>
#include <scai/benchmark/InputSetRegistry.hpp>

#include <algorithm>
#include <numeric>
#include <omp.h>

namespace bf
{

SCAI_LOG_DEF_LOGGER( Benchmark::logger, "Benchmark" );

Benchmark::Benchmark()
    : mMinTime( 0.0 ), mExecutedTime( 0 ), mActNumRep( 0 ), mNumRepitions( 1 ), mSetUpTime( 0.0 ), mTearDownTime(
        0.0 )
{
    SCAI_LOG_INFO( logger, "Benchmark()" );
}

Benchmark::Benchmark( const Benchmark& other )
    : mName( other.mName ), mGId( other.mGId ), mMinTime( other.mMinTime ), mExecutedTime( 0 ), mInputSetId(
        other.mInputSetId ), mActNumRep( 0 ), mNumRepitions( other.mNumRepitions ), mSetUpTime(
            0.0 ), mTearDownTime( 0.0 )
{
    SCAI_LOG_INFO( logger, "Benchmark( " << other.mName << ")" );
}

Benchmark::Benchmark( const std::string& id, const std::string& gid )
    : mName( id ), mGId( gid ), mMinTime( 0.0 ), mExecutedTime( 0 ), mActNumRep( 0 ), mNumRepitions( 1 ), mSetUpTime(
        0.0 ), mTearDownTime( 0.0 )
{
    SCAI_LOG_INFO( logger, "Benchmark( id = " << id << ", gid = " << gid << " )" );
}

Benchmark::Benchmark( const std::string& )
    : mMinTime( 0.0 ), mExecutedTime( 0 ), mActNumRep( 0 ), mNumRepitions( 1 ), mSetUpTime( 0.0 ), mTearDownTime(
        0.0 )
{
    SCAI_LOG_INFO( logger, "Benchmark( arg IGNORED ?)" );
}

Benchmark::Benchmark( const std::string& id, const std::string& gid, const std::string& )
    : mName( id ), mGId( gid ), mMinTime( 0.0 ), mExecutedTime( 0 ), mActNumRep( 0 ), mNumRepitions( 1 ), mSetUpTime(
        0.0 ), mTearDownTime( 0.0 )
{
    SCAI_LOG_INFO( logger, "Benchmark( id = " << id << ", gid = " << gid << ", arg = ?IGNORED? )" );
}

Benchmark::~Benchmark()
{
}

void Benchmark::run( std::ostream& out )
{
    if( doOutput() )
    {
        out << "  " << getId() << std::endl;
        out << "  " << getInputSetId() << std::endl;
    }
    initialize();
    if( doOutput() )
    {
        out << "  " << getName() << std::endl;
    }

    synchronize();
    mSetUpTime = omp_get_wtime();
    setUp();
    synchronize();
    mSetUpTime = omp_get_wtime() - mSetUpTime;

    mExecutionTimes.clear();
    double executionTime = 0.0;
    for( ; mActNumRep < mNumRepitions || mExecutedTime < mMinTime; ++mActNumRep )
    {
        synchronize();
        executionTime = omp_get_wtime();
        execute();
        synchronize();
        executionTime = omp_get_wtime() - executionTime;
        mExecutedTime += executionTime;
        mExecutionTimes.push_back( executionTime );
    }

    synchronize();
    mTearDownTime = omp_get_wtime();
    tearDown();
    synchronize();
    mTearDownTime = omp_get_wtime() - mTearDownTime;

    shutdown();
}

const std::string& Benchmark::getName() const
{
    return mName;
}

const std::string& Benchmark::getGid() const
{
    return mGId;
}

std::string Benchmark::getNumThreads() const
{
    std::stringstream threads;
    if( isThreadded() )
    {
        #pragma omp parallel
        #pragma omp master
        threads << omp_get_num_threads();
    }
    else
    {
        threads << "Unthreadded";
    }
    return threads.str();
}

bool Benchmark::isThreadded() const
{
    return false;
}

void Benchmark::setMinTime( const double minTime )
{
    mMinTime = minTime;
}

int Benchmark::getNumRepitions() const
{
    return mNumRepitions;
}

void Benchmark::setNumRepitions( const int numRepitions )
{
    mNumRepitions = numRepitions;
}

void Benchmark::setInputSetId( const std::string& inputSetId )
{
    mInputSetId = inputSetId;
}

const std::string& Benchmark::getInputSetId() const
{
    return mInputSetId;
}

double Benchmark::getSetupTime() const
{
    return mSetUpTime;
}

double Benchmark::getExecutionTime() const
{
    if( mExecutionTimes.empty() )
    {
        return 0.0;
    }
    //sort and remove the first 5% of the done measurements and the last 5% of
    //the to cut of outlier
    std::list<double> execTimes( mExecutionTimes );
    execTimes.sort();
    //cut of at least the max element and the min element if possible
    if( execTimes.size() >= 3 )
    {
        execTimes.pop_front();
        execTimes.pop_back();
    }
    while( execTimes.size() > static_cast<unsigned int>( std::max( 3, ( 9 * mNumRepitions ) / 10 ) ) )
    {
        execTimes.pop_front();
        execTimes.pop_back();
    }
    double timeSum = std::accumulate( execTimes.begin(), execTimes.end(), 0.0 );
    return timeSum / execTimes.size();
}

double Benchmark::getMinExecutionTime() const
{
    if( mExecutionTimes.empty() )
    {
        return 0.0;
    }
    return *std::min_element( mExecutionTimes.begin(), mExecutionTimes.end() );
}

double Benchmark::getMaxExecutionTime() const
{
    if( mExecutionTimes.empty() )
    {
        return 0.0;
    }
    return *std::max_element( mExecutionTimes.begin(), mExecutionTimes.end() );
}

double Benchmark::getExecutionFlops() const
{
    return getNumFloatingPointOperations() / getMinExecutionTime();
}

double Benchmark::getExecutionBandwidth() const
{
    return getProcessedBytes() / getMinExecutionTime();
}

double Benchmark::getTearDownTime() const
{
    return mTearDownTime;
}

double Benchmark::getTotalExecutionTime() const
{
    return mSetUpTime + getExecutionTime() + mTearDownTime;
}

double Benchmark::getExecutedTime() const
{
    return mExecutedTime;
}

int Benchmark::getNumActualRepititons() const
{
    return mActNumRep;
}

bool Benchmark::doOutput() const
{
    return true;
}

void Benchmark::interpreteArgs( const std::string& )
{
}

void Benchmark::synchronize() const
{
}

bool Benchmark::operator ==( const Benchmark& other )
{
    if( getId() != other.getId() )
    {
        return false;
    }
    if( getGid() != other.getGid() )
    {
        return false;
    }
    if( getName() != other.getName() )
    {
        return false;
    }
    if( getNumThreads() != other.getNumThreads() )
    {
        return false;
    }
    if( getInputSetId() != other.getInputSetId() )
    {
        return false;
    }
    return true;
}

bool Benchmark::CompareGroupIds::operator( )( const Benchmark* const arg1, const Benchmark* const arg2 )
{
    return arg1->getGid() < arg2->getGid();
}

bool Benchmark::CompareInputSetIds::operator( )( const Benchmark* const arg1, const Benchmark* const arg2 )
{
    return arg1->getInputSetId() < arg2->getInputSetId();
}

bool Benchmark::CompareValueTypeSize::operator( )( const Benchmark* const arg1, const Benchmark* const arg2 )
{
    return arg1->getValueTypeSize() < arg2->getValueTypeSize();
}

bool Benchmark::GreaterNumThreads::operator( )( const Benchmark* const arg1, const Benchmark* const arg2 )
{
    if( arg1->isThreadded() && !arg2->isThreadded() )
    {
        //arg1 is threaded and arg2 is not so num threads of arg1 are
        //considered more
        return true;
    }
    if( !arg1->isThreadded() && arg2->isThreadded() )
    {
        //arg1 is not threaded and arg2 is so num threads of arg1 are
        //considered less
        return false;
    }
    if( !arg1->isThreadded() && !arg2->isThreadded() )
    {
        //both Benchmarks are not threadded so it does not matter
        return false;
    }
    return arg1->getNumThreads() > arg2->getNumThreads();
}

bool Benchmark::LessNumThreads::operator( )( const Benchmark* const arg1, const Benchmark* const arg2 )
{
    if( arg1->isThreadded() && !arg2->isThreadded() )
    {
        //arg1 is threaded and arg2 is not so num threads of arg1 are
        //considered more
        return false;
    }
    if( !arg1->isThreadded() && arg2->isThreadded() )
    {
        //arg1 is not threaded and arg2 is so num threads of arg1 are
        //considered less
        return true;
    }
    if( !arg1->isThreadded() && !arg2->isThreadded() )
    {
        //both Benchmarks are not threadded so it does not matter
        return false;
    }
    return arg1->getNumThreads() < arg2->getNumThreads();
}

Benchmark::HasId::HasId( const std::string& id )
    : mId( id )
{
}

bool Benchmark::HasId::operator( )( const Benchmark* const benchmark )
{
    return benchmark->getId() == mId;
}

} //namespace bf
