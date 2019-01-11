/**
 * @file LAPACK/QR.cpp
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
 * @brief ToDo: Missing description in ./LAPACK/QR.cpp
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

#include "Walltime.hpp"

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
    void sgeqrf_( int* m, int* n, float* a, int* lda, float* tau, float* work,
                  int* ldwork, int* info );

    void sorgqr_( int* m, int* n, int* k, float* a,  int* lda,
                  float* tau, float* work, int* ldwork, int* info );
}

/*
    QR factorization with LAPACK routines
*/

void QR( float* Q, float* R, int m, int n )
{
    int ldwork = m * n / 4;
    int info   = 0;

    float* work = new float[ ldwork ];
    float* tau  = new float[ m ];

    sgeqrf_( &m, &n, Q, &m, tau, work, &ldwork, &info );

    printf( "SGEQRF INFO = %d\n", info );

    for ( int i = 0; i < n; ++i ) 
    {
        for ( int j = 0; j < n; ++j )
        {
            if ( j >= i )
            {
                R[index(i, j, n)] = Q[index(i, j, m)];
            }
            else
            {
                R[index(i, j, n)] = 0.0;
            }
        }
    }

    sorgqr_( &m, &n, &n, Q, &m, tau, work, &ldwork, &info );

    printf( "SORGQR INFO = %d\n", info );
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
        A[i] = 10.0f * ( float( curr ) / float( c ) - 0.5f );
    }

    return A;
}

double timestamp()
{
	return lama::Walltime::get();
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

    if ( n < 10 && m < 10 )
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

    if ( n < 10 && m < 10 )
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
                error = std::max( error, abs( A[k] - A1[k] ) );
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

