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
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <sys/time.h>

using namespace std;

// number of repetitions, for improved statistics

const int rep = 1;

#include "QR.hpp"

/* array indexing: i, j are row and column, n is the number of columns */

inline unsigned int idx( unsigned int i, unsigned j, unsigned int n )
{
    return ( j - 1 ) * n + i - 1;
}

/* fill a m x n array with random numbers, drawn by a standard
 rand32(). Used to create the reference random arrays.
 */
float* randomnumber( float* A, unsigned int seed, unsigned int m, unsigned int n )
{
    unsigned long long c = 4294967291;
    unsigned long long a = 3141592653;
    unsigned long long b = 56985;
    unsigned long long curr = seed;
    for( unsigned int i = 0; i < m * n; ++i )
    {
        curr = ( a * curr + b ) % c;
        A[i] = 10 * ( float( curr ) / c - 0.5 );
    }
    return A;
}

double timestamp()
{
    struct timeval tp;
    struct timezone tzp;

    gettimeofday( &tp, &tzp );

    return (double) tp.tv_sec + tp.tv_usec * 0.000001;
}

void printMatrix( const char* name, float* mat, unsigned int m, unsigned int n )
{
    for( int i = 1; i <= m; i++ )
    {
        printf( "%s [%d] : ", name, i );
        for( int j = 1; j <= n; j++ )
        {
            printf( " %g", mat[idx( i, j, m )] );
        }
        printf( "\n" );
    }
}

int main( int argc, char **argv )
{
    unsigned int seed, m, n;
    double start, stop;

    // read parameters from stdin
    cout << "Geben Sie den Seed ein:" << endl;
    cin >> seed;
    cout << "Geben Sie die Anzahl der Zeilen ein:" << endl;
    cin >> m;
    cout << "Geben Sie die Anzahl der Spalten ein (<= Zeilen):" << endl;
    cin >> n;
    cout << endl;

    // initialize input array

    float* A = new float[m * n];

    if( m == 3 and n == 3 )
    {
        A[idx( 1, 1, m )] = 1.0;
        A[idx( 1, 2, m )] = 2.0;
        A[idx( 1, 3, m )] = 3.0;
        A[idx( 2, 1, m )] = -1.0;
        A[idx( 2, 2, m )] = 0.0;
        A[idx( 2, 3, m )] = -3.0;
        A[idx( 3, 1, m )] = 0.0;
        A[idx( 3, 2, m )] = -2.0;
        A[idx( 3, 3, m )] = -3.0;
    }
    else
    {
        randomnumber( A, seed, m, n );
    }

    if( n < 10 and m < 10 )
    {
        printMatrix( "A", A, m, n );
    }

    // allocate result arrays for CPU

    float* Q = new float[m * n];
    float* R = new float[n * n];

    // measure calculation on the CPU

    start = timestamp();

    for( int i = 0; i < m * n; i++ )
    {
        Q[i] = A[i];
    }

    for( int i = 0; i < n * n; i++ )
    {
        R[i] = 0.0;
    }

    for( int i = 0; i < rep; ++i )
    {
        QR( Q, R, m, n );
    }

    stop = timestamp();

    cout << "timeCPU " << ( stop - start ) / rep << endl;

    double Flop = 2.0 * double( m ) * double( n ) * double( n );
    double Giga = double( 1000 ) * double( 1000 ) * double( 1000 );
    cout << "Flop = " << Flop << endl;
    double GFlops = Flop / Giga / ( stop - start );
    cout << "Performance = " << GFlops << " GFLops" << endl;

    // write out dimensions of the matrix
    cout << "N " << n << endl;
    cout << "M " << m << endl;

    if( n < 10 and m < 10 )
    {
        printMatrix( "Q", Q, m, n );
        printMatrix( "R", R, n, n );
    }

    if( n <= 512 )
    {

        // check for correct results

        float* A1 = new float[m * n];

        for( int j = 1; j <= n; j++ )
        {
            for( int i = 1; i <= m; i++ )
            {
                float S = 0.0;
                for( int k = 1; k <= j; k++ )
                {
                    S += Q[idx( i, k, m )] * R[idx( k, j, n )];
                }
                A1[idx( i, j, m )] = S;
            }
        }

        float error = 0.0;
        for( int i = 1; i <= m; i++ )
        {
            for( int j = 1; j <= n; j++ )
            {
                int k = idx( i, j, m );
                error = std::max( error, abs( A[k] - A1[k] ) );
            }
        }

        printf( "Maximal error = %g\n", error );

        // check for orthonormal base of Q

        error = 0.0;
        for( int j = 1; j <= n; j++ )
        {
            // norm(Q(:J)) should be 1.0
            float S = 0.0;
            for( int i = 1; i <= m; i++ )
            {
                S = S + Q[idx( i, j, m )] * Q[idx( i, j, m )];
            }
            error = std::max( error, abs( S - 1.0f ) );
            if( abs( S - 1.0 ) > 0.001 )
            {
                printf( "col %d of Q is not normed, val = %f\n", j, S );
            }
        }

        printf( "Norm error: %f\n", error );

        error = 0.0;

        for( int j = 1; j <= n; j++ )
        {
            for( int k = j + 1; k <= n; k++ )
            {
                // Q(:J) * Q(:K) should be 1.0
                float S = 0.0;
                for( int i = 1; i <= m; i++ )
                {
                    S = S + Q[idx( i, j, m )] * Q[idx( i, k, m )];
                }
                error = std::max( error, abs( S ) );
                if( abs( S ) > 0.001 )
                {
                    printf( "col %d %d of Q are not orthogonal, val = %f\n", j, k, S );
                }
            }
        }

        printf( "Orthogonal error: %f\n", error );

        delete[] A1;
    }

    delete[] A;
    delete[] Q;
    delete[] R;

    return 0;
}
