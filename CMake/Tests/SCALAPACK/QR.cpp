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

#include "QR.hpp"

#include <cmath>
#include <cstdio>

#ifdef HAVE_OMP
#include <omp.h>
#else
#define omp_get_num_threads() 1
#define omp_get_thread_num() 0
#endif

#include <mkl_cblas.h>

/* array indexing: i, j are row and column, n is the number of columns */

inline unsigned int index( unsigned int i, unsigned j, unsigned int n )
{
    return j * n + i;
}

/* Gram-Schmidt method on the CPU using OpenMP. Internally exchanges
 the row/column order for better performance, since most loops are
 rowwise.

 array is the m x n matrix to decompose, and the decomposed matrix
 is written to Q and R.
 */

void QR( float* Q, float* R, unsigned int m, unsigned int n )
{
    int b = 8;

    if( n > 3000 )
    {
        b = 32;
    }

    float* QS = new float[b * m];

    #pragma omp parallel
    for( unsigned int i1 = 0; i1 < n; i1 += b )
    {

        unsigned int i2 = i1 + b - 1;

        if( i2 >= n )
        {
            i2 = n - 1;
        }

        // printf("block %d-%d\n", i1, i2);

        #pragma omp single
        {
            for( unsigned i = i1; i <= i2; ++i )
            {
                cblas_sgemv( CblasColMajor, CblasTrans, m, i2 - i + 1, 1.0, &Q[index( 0, i, m )], m,
                             &Q[index( 0, i, m )], 1, 0.0, &R[index( i, i, n )], n );

                float S = 1.0 / sqrt( R[index( i, i, n )] );

                cblas_scopy( m, &Q[index( 0, i, m )], 1, &QS[index( 0, i - i1, m )], 1 );
                cblas_sscal( m, S, &Q[index( 0, i, m )], 1 );
                cblas_sscal( i2 - i + 1, S, &R[index( i, i, n )], n );
                cblas_sger( CblasColMajor, m, i2 - i, -1.0, &Q[index( 0, i, m )], 1, &R[index( i, i + 1, n )], n,
                            &Q[index( 0, i + 1, m )], m );
            }
        }

        #pragma omp for
        for( unsigned int j = i2 + 1; j < n; ++j )
        {
            for( unsigned int i = i1; i <= i2; ++i )
            {
                R[index( i, j, n )] = cblas_sdot( m, &QS[index( 0, i - i1, m )], 1, &Q[index( 0, j, m )], 1 );
                R[index( i, j, n )] /= R[index( i, i, n )];
                cblas_saxpy( m, -R[index( i, j, n )], &Q[index( 0, i, m )], 1, &Q[index( 0, j, m )], 1 );
            }
        }
    }

    delete[] QS;
}
