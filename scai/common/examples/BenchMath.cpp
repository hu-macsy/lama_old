/**
 * @file examples/BenchMath.cpp
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
 * @brief ToDo: Missing description in ./examples/BenchMath.cpp
 * @author eschricker
 * @date 12.04.2016
 */

#include <iostream>
#include <sstream>

#include <scai/common/SCAITypes.hpp>
#include <scai/common/Math.hpp>
#include <scai/common/OpenMP.hpp>
#include <scai/common/Walltime.hpp>

typedef scai::IndexType myType; // for bigger testcases need size_t

using scai::common::Math;

/*
 * Templated kernel using Math-functions from common
 */
template<typename ValueType>
ValueType asum_with_math( const myType N, const ValueType* x )
{
    ValueType asum = 0;
    #pragma omp parallel
    {
        ValueType thread_asum = 0;
        #pragma omp for

        for ( myType i = 0; i < N; ++i )
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
double asum_explicit( const myType N, const double* x )
{
    double asum = 0;
    #pragma omp parallel for reduction(+:asum)

    for ( myType i = 0; i < N; ++i )
    {
        asum += std::abs( x[i] );
    }

    return asum;
}

/*
 * Function for setting array to a specific value
 */
template<typename ValueType>
void set( const myType N, ValueType* x, const ValueType val );

/*
 * Calculate average time
 */
double avg_time( const myType N, const double* time_data );

/*
 * Parse commandline arguments
 */
void parseArguments( int argc, char** argv, myType& repitions, myType& N );

/*
 * Main
 */
int main( int argc, char** argv )
{
    myType repitions, N;
    /*
     * Default values
     */
    repitions = 30;
    N = 1 << 28; // 2^28 ~ 2GB of data

    if ( argc > 1 )
    {
        parseArguments( argc, argv, repitions, N );
    }

    std::cout << "Repitions: " << repitions << std::endl;
    std::cout << "Number of elements: " << N << std::endl;
    /*
     * Creation of data array
     */
    double* d = new double[N];
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
    for ( myType i = 0; i < repitions; ++i )
    {
        tmpTime = scai::common::Walltime::get();
        asum_with_math( N, d );
        mathTime[i] = scai::common::Walltime::get() - tmpTime;
    }

    std::cout << "Math avg: " << avg_time( repitions, mathTime ) << std::endl;

    /*
     * Explicit kernel
     */
    for ( myType i = 0; i < repitions; ++i )
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
void set( const myType N, ValueType* x, const ValueType val )
{
    #pragma omp parallel for

    for ( myType i = 0; i < N; ++i )
    {
        x[i] = val;
    }
}

double avg_time( const myType N, const double* time_data )
{
    double min, max, avg;
    min = max = time_data[1];
    avg = 0;

    for ( myType i = 1; i < N - 1; ++i )
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

void parseArguments( int argc, char** argv, myType& repitions, myType& N )
{
    for ( int i = 1; i < argc; ++i )
    {
        std::stringstream s( argv[i] );
        std::string content;

        while ( s >> content )
        {
            if ( content == "-n" || content == "--nr-elements" )
            {
                // set number of elemens
                s.clear();
                s.str( argv[++i] );
                s >> N;
            }
            else if ( content == "-r" || content == "--repitions" )
            {
                // set repitions
                s.clear();
                s.str( argv[++i] );
                s >> repitions;

                if ( repitions < 5 )
                {
                    std::cout << "Minimum repitions is 5" << std::endl;
                }
            }
            else if ( content == "-h" || content == "--help" )
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
