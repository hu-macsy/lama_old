/**
 * @file directSolver.cpp
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
 * @brief Demo of using Intel PARDISO by LAMA
 * @author Thomas Brandes
 * @date 01.07.2016
 */

#include <scai/lama.hpp>
#include <scai/hmemo.hpp>

#include <mkl.h>

using namespace scai;
using namespace lama;

/**
 *  Main program
 *
 *  - first arg is filename for input matrix
 *  - all other arguments are passed to the configuration lamaconf
 *  - configuration will contain all information to setup the solver for the input matrix
 */
int main( int argc, const char* argv[] )
{
    if ( argc != 2 )
    {
        return -1;
    }

    std::string matrixFilename = argv[1];

    CSRSparseMatrix<double> matrix;

    matrix.readFromFile( matrixFilename );

    std::cout << "Matrix = " << matrix << std::endl;

    DenseVector<double> sol( matrix.getNumColumns(), 1.0 );

    DenseVector<double> rhs( matrix * sol );

    // now call Pardiso

    CSRStorage<double>& csr = const_cast<CSRStorage<double>&>( matrix.getLocalStorage() );

    hmemo::WriteAccess<IndexType> ia( csr.getIA() );
    hmemo::WriteAccess<IndexType> ja( csr.getJA() );
    hmemo::WriteAccess<double> a( csr.getValues() );
    hmemo::WriteAccess<double> b( rhs.getLocalValues() );
    hmemo::WriteAccess<double> x( sol.getLocalValues() );

    MKL_INT n = csr.getNumRows();

    printf( "n = %d\n", n );

    /*

    MKL_INT n = 5;
    MKL_INT ia[ 6] = { 0, 3, 5, 8, 11, 13 };
    MKL_INT ja[13] = { 0, 1, 3,
        0, 1,
        2, 3, 4,
        0, 2, 3,
        1, 4 };
    double a[18] = { 1.0, -1.0, -3.0,
        -2.0, 5.0,
        4.0, 6.0, 4.0,
        -4.0, 2.0, 7.0,
        8.0, -5.0 };

    // RHS and solution vectors.

    double b[5], x[5];

    */

    MKL_INT mtype = 11; // Real unsymmetric matrix
    MKL_INT nrhs  = 1;  // Number of right hand sides

    /* Internal solver memory pointer pt, */
    /* 32-bit: int pt[64]; 64-bit: long int pt[64] */
    /* or void *pt[64] should be OK on both architectures */

    void* pt[64];

    MKL_INT iparm[64];  // control parameters

    MKL_INT maxfct, mnum, phase, error, msglvl;
    /* Auxiliary variables. */
    MKL_INT i;
    double ddum; /* Double dummy */
    MKL_INT idum; /* Integer dummy. */

    /* -------------------------------------------------------------------- */
    /* .. Setup Pardiso control parameters. */
    /* -------------------------------------------------------------------- */
    for ( i = 0; i < 64; i++ )
    {
        iparm[i] = 0;
    }

    iparm[0] = 1; /* No solver default */
    iparm[1] = 2; /* Fill-in reordering from METIS */
    /* Numbers of processors, value of OMP_NUM_THREADS */
    iparm[2] = 1;
    iparm[3] = 0; /* No iterative-direct algorithm */
    iparm[4] = 0; /* No user fill-in reducing permutation */
    iparm[5] = 0; /* Write solution into x */
    iparm[6] = 0; /* Not in use */
    iparm[7] = 2; /* Max numbers of iterative refinement steps */
    iparm[8] = 0; /* Not in use */
    iparm[9] = 13; /* Perturb the pivot elements with 1E-13 */
    iparm[10] = 1; /* Use nonsymmetric permutation and scaling MPS */
    iparm[11] = 0; /* Not in use */
    iparm[12] = 1; /* Maximum weighted matching algorithm is switched-on (default for non-symmetric) */
    iparm[13] = 0; /* Output: Number of perturbed pivots */
    iparm[14] = 0; /* Not in use */
    iparm[15] = 0; /* Not in use */
    iparm[16] = 0; /* Not in use */
    iparm[17] = -1; /* Output: Number of nonzeros in the factor LU */
    iparm[18] = -1; /* Output: Mflops for LU factorization */
    iparm[19] = 0; /* Output: Numbers of CG Iterations */
    iparm[34] = 1; /* PARDISO use C-style indexing for ia and ja arrays */
    maxfct = 1; /* Maximum number of numerical factorizations. */
    mnum = 1; /* Which factorization to use. */
    msglvl = 1; /* Print statistical information in file */
    error = 0; /* Initialize error flag */

    /* -------------------------------------------------------------------- */
    /* .. Initialize the internal solver memory pointer. This is only */
    /* necessary for the FIRST call of the PARDISO solver. */
    /* -------------------------------------------------------------------- */
    for ( i = 0; i < 64; i++ )
    {
        pt[i] = 0;
    }

    /* -------------------------------------------------------------------- */
    /* .. Reordering and Symbolic Factorization. This step also allocates */
    /* all memory that is necessary for the factorization. */
    /* -------------------------------------------------------------------- */
    phase = 11;

    printf( "call PARDISO phase11\n" );

    PARDISO ( pt, &maxfct, &mnum, &mtype, &phase,
              &n, a, ia, ja, &idum, &nrhs,
              iparm, &msglvl, &ddum, &ddum, &error );

    if ( error != 0 )
    {
        printf( "\nERROR during symbolic factorization: %d", error );
        exit( 1 );
    }

    printf( "\nReordering completed ... " );
    printf( "\nNumber of nonzeros in factors = %d", iparm[17] );
    /* -------------------------------------------------------------------- */
    /* .. Numerical factorization. */
    /* -------------------------------------------------------------------- */
    phase = 22;
    PARDISO ( pt, &maxfct, &mnum, &mtype, &phase,
              &n, a, ia, ja, &idum, &nrhs,
              iparm, &msglvl, &ddum, &ddum, &error );

    if ( error != 0 )
    {
        printf( "\nERROR during numerical factorization: %d", error );
        exit( 2 );
    }

    printf( "\nFactorization completed ... " );
    /* -------------------------------------------------------------------- */
    /* .. Back substitution and iterative refinement. */
    /* -------------------------------------------------------------------- */
    phase = 33;
    iparm[7] = 2; /* Max numbers of iterative refinement steps. */

    /* Set right hand side to one. */
    for ( i = 0; i < n; i++ )
    {
        b[i] = 1;
    }

    PARDISO ( pt, &maxfct, &mnum, &mtype, &phase,
              &n, a, ia, ja, &idum, &nrhs,
              iparm, &msglvl, b, x, &error );

    if ( error != 0 )
    {
        printf( "\nERROR during solution: %d", error );
        exit( 3 );
    }

    printf( "\nSolve completed ... " );
    printf( "\nThe solution of the system is: " );

    for ( i = 0; i < n; i++ )
    {
        printf( "\n x [%d] = % f", i, x[i] );
    }

    printf ( "\n" );
    /* -------------------------------------------------------------------- */
    /* .. Termination and release of memory. */
    /* -------------------------------------------------------------------- */
    phase = -1; /* Release internal memory. */
    PARDISO ( pt, &maxfct, &mnum, &mtype, &phase,
              &n, &ddum, ia, ja, &idum, &nrhs,
              iparm, &msglvl, &ddum, &ddum, &error );
    return 0;

}
