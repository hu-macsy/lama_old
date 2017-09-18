/**
 * @file Benchmark.cpp
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
 * @brief Implementation of methods for class Benchmark.
 * @author jiri
 * @date 04.05.2010
 */

#include <scai/benchmark/Benchmark.hpp>
#include <scai/benchmark/Parser.hpp>

#include <scai/common/Walltime.hpp>
#include <scai/common/OpenMP.hpp>

#include <algorithm>
#include <numeric>

namespace scai
{

namespace benchmark
{

SCAI_LOG_DEF_LOGGER( Benchmark::logger, "Benchmark" );

Benchmark::Benchmark( const std::string& name, const std::string& group ) : 

    mName( name ), 
    mGroup( group ), 
    mMinTime( 0.0 ), 
    mExecutedTime( 0 ), 
    mActNumRep( 0 ), 
    mNumRepitions( 1 ), 
    mSetUpTime( 0.0 ), 
    mTearDownTime( 0.0 )
{
    SCAI_LOG_INFO( logger, "Benchmark( name = " << name << ", group = " << group << " )" );
}

Benchmark::~Benchmark()
{
}

void Benchmark::writeAt( std::ostream& stream ) const
{
    stream << "Benchmark( name = " << mName << ", group = " << mGroup << " )";
}

void Benchmark::run( std::ostream& out )
{
    if ( doOutput() )
    {
        out << "  " << getCreateId() << std::endl;
        out << "  " << getInputSetId() << std::endl;
    }

    initialize();

    if ( doOutput() )
    {
        out << "  " << getName() << std::endl;
    }

    synchronize();
    mSetUpTime = common::Walltime::get();
    setUp();
    synchronize();
    mSetUpTime = common::Walltime::get() - mSetUpTime;

    mExecutionTimes.clear();
    double executionTime = 0.0;

    for ( ; mActNumRep < mNumRepitions || mExecutedTime < mMinTime; ++mActNumRep )
    {
        synchronize();
        executionTime = common::Walltime::get();
        execute();
        synchronize();
        executionTime = common::Walltime::get() - executionTime;
        mExecutedTime += executionTime;
        mExecutionTimes.push_back( executionTime );
    }

    synchronize();
    mTearDownTime = common::Walltime::get();
    tearDown();
    synchronize();
    mTearDownTime = common::Walltime::get() - mTearDownTime;

    shutdown();
}

const std::string& Benchmark::getName() const
{
    return mName;
}

const std::string& Benchmark::getGroup() const
{
    return mGroup;
}

std::string Benchmark::getNumThreads() const
{
    std::stringstream threads;

    if ( isThreadded() )
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
    SCAI_ASSERT_GT_ERROR( numRepitions, 0, "#repitions must be positive" )

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
    if ( mExecutionTimes.empty() )
    {
        return 0.0;
    }

    //sort and remove the first 5% of the done measurements and the last 5% of
    //the to cut of outlier
    std::list<double> execTimes( mExecutionTimes );
    execTimes.sort();

    //cut of at least the max element and the min element if possible
    if ( execTimes.size() >= 3 )
    {
        execTimes.pop_front();
        execTimes.pop_back();
    }

    while ( execTimes.size() > static_cast<unsigned int>( std::max( 3, ( 9 * mNumRepitions ) / 10 ) ) )
    {
        execTimes.pop_front();
        execTimes.pop_back();
    }

    double timeSum = std::accumulate( execTimes.begin(), execTimes.end(), 0.0 );
    return timeSum / execTimes.size();
}

double Benchmark::getMinExecutionTime() const
{
    if ( mExecutionTimes.empty() )
    {
        return 0.0;
    }

    return *std::min_element( mExecutionTimes.begin(), mExecutionTimes.end() );
}

double Benchmark::getMaxExecutionTime() const
{
    if ( mExecutionTimes.empty() )
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

bool Benchmark::CompareGroupIds::operator( )( const Benchmark* const arg1, const Benchmark* const arg2 )
{
    return arg1->getGroup() < arg2->getGroup();
}

bool Benchmark::CompareInputSetIds::operator( )( const Benchmark* const arg1, const Benchmark* const arg2 )
{
    return arg1->getInputSetId() < arg2->getInputSetId();
}

bool Benchmark::CompareValueTypeSize::operator( )( const Benchmark* const arg1, const Benchmark* const arg2 )
{
    return   common::typeSize( arg1->getValueType() )
           < common::typeSize( arg2->getValueType() );
}

bool Benchmark::GreaterNumThreads::operator( )( const Benchmark* const arg1, const Benchmark* const arg2 )
{
    if ( arg1->isThreadded() && !arg2->isThreadded() )
    {
        //arg1 is threaded and arg2 is not so num threads of arg1 are
        //considered more
        return true;
    }

    if ( !arg1->isThreadded() && arg2->isThreadded() )
    {
        //arg1 is not threaded and arg2 is so num threads of arg1 are
        //considered less
        return false;
    }

    if ( !arg1->isThreadded() && !arg2->isThreadded() )
    {
        //both Benchmarks are not threadded so it does not matter
        return false;
    }

    return arg1->getNumThreads() > arg2->getNumThreads();
}

bool Benchmark::LessNumThreads::operator( )( const Benchmark* const arg1, const Benchmark* const arg2 )
{
    if ( arg1->isThreadded() && !arg2->isThreadded() )
    {
        //arg1 is threaded and arg2 is not so num threads of arg1 are
        //considered more
        return false;
    }

    if ( !arg1->isThreadded() && arg2->isThreadded() )
    {
        //arg1 is not threaded and arg2 is so num threads of arg1 are
        //considered less
        return true;
    }

    if ( !arg1->isThreadded() && !arg2->isThreadded() )
    {
        //both Benchmarks are not threadded so it does not matter
        return false;
    }

    return arg1->getNumThreads() < arg2->getNumThreads();
}

Benchmark* Benchmark::createWithArgument( const std::string& specification )
{
    std::string keyValue;
    std::string argument;

    parseCommand( keyValue, argument, specification );

    return create( keyValue, argument );
}

} //namespace benchmark

} // scai
