/*
 * BenchMath.cpp
 *
 *  Created on: Apr 12, 2016
 *      Author: eschricker
 */

#include <iostream>

#include <scai/common/SCAITypes.hpp>
#include <scai/common/Complex.hpp>
#include <scai/common/Math.hpp>
#include <scai/common/OpenMP.hpp>
#include <scai/common/Walltime.hpp>

using scai::common::Math;

template<typename ValueType>
ValueType asum_with_math( const size_t N, const ValueType* x )
{
    ValueType asum = 0;
    #pragma omp parallel
    {
        ValueType thread_asum = 0;

        #pragma omp for
        for( size_t i = 0; i < N; ++i)
        {
            thread_asum += Math::abs( Math::real( x[i] ) + Math::imag( x[i] ) );
        }

        atomicAdd( asum, thread_asum );
    }

    return asum;
}

double asum_explicit( const size_t N, const double* x )
{
    double asum = 0;

    #pragma omp parallel for reduction(+:asum)
    for( size_t i = 0; i < N; ++i)
    {
        asum += std::abs( x[i] );
    }

    return asum;
}

template<typename ValueType>
void set( const size_t N, ValueType* x, const ValueType val )
{
#pragma omp parallel for
    for( size_t i = 0; i < N; ++i )
    {
        x[i] = val;
    }
}

double avg_time( const size_t N, const double* time_data )
{
    double min, max, avg;
    min = max = time_data[1];
    avg = 0;
    for( size_t i = 1; i < N-1; ++i )
    {
        min = Math::min( min, time_data[i] );
        max = Math::max( max, time_data[i] );
        avg += time_data[i];
    }

    avg -= min;
    avg -= max;
    avg /= ( N - 4 );

    return avg;
}

int main()
{
    const size_t repitions = 30;
    const size_t N = 1 << 28; // 2^28 ~ 2GB of data
    double *d = new double[N];

    double tmpTime;
    double* mathTime = new double[repitions];
    double* explicitTime = new double[repitions];

    set( N, d, 3.0 );

    for( size_t i = 0; i < repitions; ++i )
    {
        tmpTime = scai::common::Walltime::get();
        asum_with_math( N, d );
        mathTime[i] = scai::common::Walltime::get() - tmpTime;
    }
    std::cout << "Math avg: " << avg_time( repitions, mathTime ) << std::endl;

    for( size_t i = 0; i < repitions; ++i )
    {
        tmpTime = scai::common::Walltime::get();
        asum_explicit( N, d );
        explicitTime[i] = scai::common::Walltime::get() - tmpTime;
    }
    std::cout << "Explicit avg: " << avg_time( repitions, explicitTime ) << std::endl;

    delete[] d;
    delete[] mathTime;
    delete[] explicitTime;
}
