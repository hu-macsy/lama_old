/*
 * BenchMath.cpp
 *
 *  Created on: Apr 12, 2016
 *      Author: eschricker
 */

#include <iostream>
#include <sstream>

#include <scai/common/SCAITypes.hpp>
#include <scai/common/Complex.hpp>
#include <scai/common/Math.hpp>
#include <scai/common/OpenMP.hpp>
#include <scai/common/Walltime.hpp>

using scai::common::Math;

/*
 * Templated kernel using Math-functions from common
 */
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

/*
 * Explicit kernel for double
 */
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

/*
 * Function for setting array to a specific value
 */
template<typename ValueType>
void set( const size_t N, ValueType* x, const ValueType val );

/*
 * Calculate average time
 */
double avg_time( const size_t N, const double* time_data );

/*
 * Parse commandline arguments
 */
void parseArguments( int argc, char **argv, size_t& repitions, size_t& N );

/*
 * Main
 */
int main(int argc, char **argv)
{
    size_t repitions, N;

    /*
     * Default values
     */
    repitions = 30;
    N = 1 << 28; // 2^28 ~ 2GB of data

    if( argc > 1 )
    {
        parseArguments( argc, argv, repitions, N );
    }

    std::cout << "Repitions: " << repitions << std::endl;
    std::cout << "Number of elements: " << N << std::endl;

    /*
     * Creation of data array
     */
    double *d = new double[N];

    /*
     * Time measurement
     */
    double tmpTime;
    double* mathTime = new double[repitions];
    double* explicitTime = new double[repitions];

    set( N, d, 3.0 );

    /*
     * Math kernel
     */
    for( size_t i = 0; i < repitions; ++i )
    {
        tmpTime = scai::common::Walltime::get();
        asum_with_math( N, d );
        mathTime[i] = scai::common::Walltime::get() - tmpTime;
    }
    std::cout << "Math avg: " << avg_time( repitions, mathTime ) << std::endl;

    /*
     * Explicit kernel
     */
    for( size_t i = 0; i < repitions; ++i )
    {
        tmpTime = scai::common::Walltime::get();
        asum_explicit( N, d );
        explicitTime[i] = scai::common::Walltime::get() - tmpTime;
    }
    std::cout << "Explicit avg: " << avg_time( repitions, explicitTime ) << std::endl;

    /*
     * Memory clean up
     */
    delete[] d;
    delete[] mathTime;
    delete[] explicitTime;
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

void parseArguments( int argc, char **argv, size_t& repitions, size_t& N )
{
    for( int i = 0; i < argc; ++i )
    {
        std::stringstream s( argv[i] );
        std::string content;

        while( s >> content )
        {
            if( content == "-n" || content == "--nr-elements" )
            {
                // set number of elemens

                s.clear();
                s.str( argv[++i] );
                s >> N;
            }
            else if( content == "-r" || content == "--repitions" )
            {
                // set repitions
                s.clear();
                s.str( argv[++i] );
                s >> repitions;

                if( repitions < 5 )
                {
                    std::cout << "Minimum repitions is 5" << std::endl;
                }
            }
            else if( content == "-h" || content == "--help" )
            {
                // Display help
                std::cout << argv[0] << " [options]" << std::endl << std::endl <<
                                "Valid options:" << std::endl << std::endl <<
                                "\t-h --help\t\t show this help" << std::endl <<
                                "\t-r --repitions\t\t specify number of repitions" << std::endl <<
                                "\t-n --nr-elements\t specify number of elements" << std::endl;
                exit( 0 );
            }
            else
            {
                std::cout << "Unknown value " << content << std::endl;
            }
        }
    }
}