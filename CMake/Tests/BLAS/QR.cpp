/**
 * @file BLAS/QR.cpp
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
 * @brief Example QR decomposition.
 * @author Jan Ecker
 * @date 20.03.2013
 */

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
#include <algorithm>
#include <sys/time.h>

using namespace std;

// number of repetitions, for improved statistics

const int rep = 1;

/* array indexing: i, j are row and column, n is the number of columns */

inline int index( int i, int j, int n )
{
    return j * n + i;
}

extern "C"
{
    void sgemv_( char* transa, int* m, int* n, float* alpha, float* a, int* lda,
                 float* x, int* incx, float* beta, float* y, int* incy );

    void sscal_( int* n, float* alpha, float* x, int* incx );

    void sger_( int* m, int* n, float* alpha, float* x, int* incx,
                float* y, int* incy, float* a, int* lda );
}

/* Gram-Schmidt method on the CPU using OpenMP. Internally exchanges
   the row/column order for better performance, since most loops are
   rowwise.

   array is the m x n matrix to decompose, and the decomposed matrix
   is written to Q and R.
*/

void QR( float* Q, float* R, int m, int n )
{
    for ( int k = 0; k < n; ++k )
    {
        char kind = 'T';
        float f_one = 1.0;
        float f_mone = -1.0;
        float f_zero = 0.0;
        int one = 1;
        int n1 = n - k;
        sgemv_( &kind, &m, &n1, &f_one, Q + index( 0, k, m ), &m, Q + index( 0, k, m ), &one, &f_zero, R + index( k, k, n ), &n );

        float S = sqrt( R[index( k, k, n )] );
        float invS = 1.0 / S;

        sscal_( &m, &invS, Q + index( 0, k, m ), &one);

        sscal_( &n1, &invS, R + index( k, k, n), &n );

        // call sger(M, N-K, -1.0, Q(1,K), 1, R(K,K+1), N, Q(1,K+1), M)

        n1 = n - k - 1;
        sger_( &m, &n1, &f_mone, Q + index(0, k, m), &one, R + index(k, k+1, n ), &n,
                                 Q + index(0, k+1, m), &m );
    }
}

/* array indexing: i, j are row and column, n is the number of columns */

inline int idx( int i, int j, int n )
{
    return ( j - 1 ) * n + i - 1;
}
/* fill a m x n array with random numbers, drawn by a standard
   rand32(). Used to create the reference random arrays.
 */
float* randomnumber( float* A, int seed, int m, int n )
{
    unsigned long long c = 4294967291;
    unsigned long long a = 3141592653;
    unsigned long long b = 56985;
    unsigned long long curr = seed;

    for ( int i = 0; i < m * n; ++i )
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
    return( double ) tp.tv_sec + tp.tv_usec * 0.000001;
}

void printMatrix( const char* name, float* mat, int m, int n )
{
    for ( int i = 1; i <= m; i++ )
    {
        printf( "%s [%d] : ", name, i );

        for ( int j = 1; j <= n; j++ )
        {
            printf( " %g", mat[idx( i, j, m )] );
        }

        printf ( "\n" );
    }
}

int main( int argc, char** argv )
{
    int seed = 11591;
    int m;
    int n;
    double start, stop;
    // read parameters from stdin
    cout << "Geben Sie die Anzahl der Zeilen ein:" << endl;
    cin >> m;
    cout << "Geben Sie die Anzahl der Spalten ein (<= Zeilen):" << endl;
    cin >> n;
    cout << endl;
    // initialize input array
    float* A = new float[m * n];

    randomnumber( A, seed, m, n );

    if ( n < 10 and m < 10 )
    {
        printMatrix( "A", A, m, n );
    }

    // allocate result arrays for CPU
    float* Q = new float [m * n];
    float* R = new float [n * n];
    // measure calculation on the CPU
    start = timestamp();

    for ( int i = 0; i < m * n; i++ )
    {
        Q[i] = A[i];
    }

    for ( int i = 0; i < n * n; i++ )
    {
        R[i] = 0.0;
    }

    for ( int i = 0; i < rep; ++i )
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

    if ( n < 10 and m < 10 )
    {
        printMatrix( "Q", Q, m, n );
        printMatrix( "R", R, n, n );
    }

    if ( n <= 512 )
    {
        // check for correct results
        float* A1 = new float[m * n];

        for ( int j = 1; j <= n; j++ )
        {
            for ( int i = 1; i <= m; i++ )
            {
                float S = 0.0;

                for ( int k = 1; k <= j; k++ )
                {
                    S += Q[idx( i, k, m )] * R[idx( k, j, n )];
                }

                A1[idx( i, j, m )] = S;
            }
        }

        float error = 0.0;

        for ( int i = 1; i <= m; i++ )
        {
            for ( int j = 1; j <= n; j++ )
            {
                int k = idx( i, j, m );
                float diff = A[k] - A1[k];
                diff = abs( diff );
                error = std::max( error, diff );
            }
        }

        printf( "Maximal error = %g\n", error );
        delete [] A1;
    }

    delete [] A;
    delete [] Q;
    delete [] R;
    return 0;
}

