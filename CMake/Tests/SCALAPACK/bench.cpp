/*
 Copyright (C) 2010 Institute for Computational Physics, University of Stuttgart, Pfaffenwaldring 27, 70569 Germany
 Copyright (C) 2009 Fraunhofer SCAI, Schloss Birlinghoven, 53754 Sankt Augustin, Germany

 This file is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 This file is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <sys/time.h>

using namespace std;

// number of repetitions, for improved statistics

const int rep = 1;

#include "QR.hpp"

/* fill a m x n array with random numbers, drawn by a standard
 rand32(). Used to create the reference random arrays.
 */

void randomnumber( float* array, unsigned int seed, unsigned int m, unsigned int n )
{
    unsigned long long c = 42949;
    unsigned long long a = 31415;
    unsigned long long b = 56985;
    unsigned long long curr = seed;
    for( unsigned int i = 0; i < m * n; ++i )
    {
        curr = ( a * curr + b ) % c;
        // printf("curr = %lld\n", curr);
        array[i] = 10 * ( float( curr ) / c - 0.5 );
    }
}

double timestamp()
{
    struct timeval tp;
    struct timezone tzp;

    gettimeofday( &tp, &tzp );

    return (double) tp.tv_sec + tp.tv_usec * 0.000001;
}

#define IDX2F(i,j,ld) ((((j)-1)*(ld))+((i)-1))

int main( int argc, char **argv )
{
    unsigned int seed, m, n;

    int sizes[] =
    { 512, 1024, 2048, 3072, 4096, 8192 };

    for( int i = 0; i < sizeof( sizes ) / sizeof(int); i++ )
    {

        double time = 0; // sums up timing for each run

        seed = 1514;
        m = sizes[i];
        n = sizes[i];

        for( int k = 0; k < rep; k++ )
        {

            float* A = new float[m * n];

            randomnumber( A, seed, m, n );

            // allocate result arrays for GPU

            float* Q = new float[m * n];
            float* R = new float[n * n];

            for( int i = 0; i < n * n; i++ )
            {
                R[i] = 0.0;
            }
            for( int i = 0; i < m * n; i++ )
            {
                Q[i] = A[i];
            }

            double start = timestamp();

            QR( Q, R, m, n );

            double time1 = timestamp() - start;
            ;

            cout << "time run " << k + 1 << ": " << time1 << " seconds" << endl;

            time += time1;

            delete[] A;
            delete[] Q;
            delete[] R;
        }

        cout << "time " << time << endl;

        // write out dimensions of the matrix

        cout << "M = " << m << ", N = " << n << endl;

        double Flop = double( n ) * double( n ) * double( m ) * 2.0 * rep;
        double GFlops = Flop * 0.001 * 0.001 * 0.001;
        GFlops = GFlops / time;

        cout << "GFlops = " << GFlops << endl;

    }

    return 0;
}
