/**
 * @file OpenMPSCALAPACK.cpp
 *
 * @license
 * Copyright (c) 2012
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 * @endlicense
 *
 * @brief OpenMPSCALAPACK.cpp
 * @author lschubert
 * @date 04.07.2012
 * $Id$
 */

// hpp
#include <lama/openmp/OpenMPSCALAPACK.hpp>

// others
#include <lama/openmp/SCALAPACKHelper.hpp>

#include <lama/Communicator.hpp>
#include <lama/LAMAInterface.hpp>
#include <lama/LAMAInterfaceRegistry.hpp>

// macros
#include <lama/macros/unused.hpp>

// boost
#include <boost/scoped_array.hpp>

#include <mkl_scalapack.h>
#include <mkl_blacs.h>

#ifndef WIN32
#define F77_descinit        descinit_
#define F77_numroc          numroc_
#define F77_blacs_get       blacs_get_
#define F77_blacs_gridinit  blacs_gridinit_
#define F77_blacs_gridinfo  blacs_gridinfo_
#define F77_blacs_gridexit  blacs_gridexit_
#define F77_blacs_freebuff  blacs_freebuff_
#else
#define F77_descinit        DESCINIT
#define F77_numroc          NUMROC
#define F77_blacs_get       BLACS_GET
#define F77_blacs_gridinit  BLACS_GRIDINIT
#define F77_blacs_gridinfo  BLACS_GRIDINFO
#define F77_blacs_gridexit  BLACS_GRIDEXIT
#define F77_blacs_freebuff  BLACS_FREEBUFF
#endif

using namespace lama;

namespace lama
{

LAMA_LOG_DEF_LOGGER( OpenMPSCALAPACK::logger, "OpenMP.SCALAPACK" )

template<>
IndexType OpenMPSCALAPACK::pgetrf(
    const IndexType m,
    const IndexType n,
    const float* const A,
    const IndexType ia,
    const IndexType ja,
    IndexType* const descA,
    IndexType* const ipiv )
{
    IndexType info = -1;
    psgetrf( (int*) &m, (int*) &n, (float*) A, (int*) &ia, (int*) &ja, (int*) descA, ipiv, &info );
    return info;
}

template<>
IndexType OpenMPSCALAPACK::pgetrf(
    const IndexType m,
    const IndexType n,
    const double* const A,
    const IndexType ia,
    const IndexType ja,
    IndexType* const descA,
    IndexType* const ipiv )
{
    IndexType info = -1;
    pdgetrf( (int*) &m, (int*) &n, (double*) A, (int*) &ia, (int*) &ja, (int*) descA, ipiv, &info );
    return info;
}

template<>
IndexType OpenMPSCALAPACK::pgetri(
    const IndexType n,
    const float* const A,
    const IndexType ia,
    const IndexType ja,
    IndexType* const descA,
    IndexType* const ipiv,
    const float* const work,
    IndexType lwork,
    IndexType* const iwork,
    IndexType liwork )
{
    IndexType info = -1;
    // call of Fortran/C routine, needs always pointers
    psgetri( (int*) &n, (float*) A, (int*) &ia, (int*) &ja, (int*) descA, ipiv, (float*) work, (int*) &lwork, iwork,
             (int*) &liwork, &info );
    return info;
}

template<>
IndexType OpenMPSCALAPACK::pgetri(
    const IndexType n,
    const double* const A,
    const IndexType ia,
    const IndexType ja,
    IndexType* const descA,
    IndexType* const ipiv,
    const double* const work,
    IndexType lwork,
    IndexType* const iwork,
    IndexType liwork )
{
    IndexType info = -1;
    // call of Fortran/C routine, needs always pointers
    pdgetri( (int*) &n, (double*) A, (int*) &ia, (int*) &ja, (int*) descA, ipiv, (double*) work, (int*) &lwork, iwork,
             (int*) &liwork, &info );
    return info;
}

template<typename T>
void OpenMPSCALAPACK::inverse( const IndexType n, const IndexType nB, const T* a, const class Communicator& comm )
{
    int contxt = -1;

    SCALAPACKHelper::BLACS_Context( &contxt );

    LAMA_LOG_INFO( logger, "Default system context = " << contxt )

    // Some tricky stuff here: we compute inverse of the transposed matrix
    // So we do as if the matrix is distributed among the columns

    const int m = n; // square matrix

    const int npRow = 1;
    const int npCol = comm.getSize();

    //After this point assume that we invert a matrix with m rows an n columns
    //that is stored column major order, so no further conversions are necessary
    //to call fortran routines

    int myRow;
    int myCol;

    SCALAPACKHelper::GridInit( &contxt, "Col-major", npRow, npCol );
    SCALAPACKHelper::GridInfo( contxt, npRow, npCol, &myRow, &myCol );

    LAMA_LOG_INFO( logger, comm << ": me is " << myRow << " x " << myCol << " of " << npRow << " x " << npCol )

    // Note: numbering must be the same as communicator

    LAMA_ASSERT_EQUAL_ERROR( myCol, comm.getRank() )

    int descA[9]; // descriptor needed for matrix A

    const int mb_a = nB; // block size in first dimension
    const int nb_a = mb_a; // block size in second dimension

    const int irsrc = 0; // processor with first row
    const int icsrc = 0; // processor with first col
    const int lld = m; // max (1, np) + imidpapd

    int ierr = SCALAPACKHelper::DescInit( descA, m, n, mb_a, nb_a, irsrc, icsrc, contxt, lld );

    LAMA_LOG_INFO( logger, "lama_DescInit, error = " << ierr )

    if ( ierr != 0 )
    {
        LAMA_THROWEXCEPTION( "lama_DescInit failed, error = " << ierr )
    }

    const int mp = SCALAPACKHelper::numRoC( m, mb_a, myRow, npRow );

    LAMA_LOG_INFO( logger, "mp = " << mp << ", np = " << SCALAPACKHelper::numRoC( n, nb_a, myCol, npCol ) )

    boost::scoped_array<IndexType> permutation( new IndexType[mp + mb_a] );

    LAMA_LOG_INFO( logger, "now call pgetrf" )

    int info = pgetrf( m, n, a, 1, 1, descA, permutation.get() );

    LAMA_LOG_INFO( logger, "pgetrf: info = " << info )

    if ( info != 0 )
    {
        LAMA_THROWEXCEPTION( "pgetrf failed, error = " << info )
    }

    int lwork = -1;
    int liwork = -1;

    T workTmp = 0;
    int iworkTmp = 0;

    //Let scalapack compute the optimal values of lwork and liwork

    info = pgetri( m, a, 1, 1, descA, permutation.get(), &workTmp, lwork, &iworkTmp, liwork );

    if ( info != 0 )
    {
        LAMA_THROWEXCEPTION( "pgetri (first call) failed, error = " << info )
    }

    lwork = static_cast<int>( workTmp );
    liwork = iworkTmp;

    boost::scoped_array<T> work( new T[lwork] );
    boost::scoped_array<int> iwork( new int[liwork] );

    info = pgetri( m, a, 1, 1, descA, permutation.get(), work.get(), lwork, iwork.get(), liwork );

    if ( info != 0 )
    {
        LAMA_THROWEXCEPTION( "pgetri (final call) failed, error = " << info )
    }
}

/* --------------------------------------------------------------------------- */
/*     Template instantiations via registration routine                        */
/* --------------------------------------------------------------------------- */

void OpenMPSCALAPACK::setInterface( BLASInterface& BLAS )
{
    // Note: macro takes advantage of same name for routines and type definitions 
    //       ( e.g. routine CUDABLAS1::sum<T> is set for BLAS::BLAS1::sum variable

    LAMA_INTERFACE_REGISTER_T( BLAS, inverse, float )
    LAMA_INTERFACE_REGISTER_T( BLAS, inverse, double )

    // other routines are not used by LAMA yet
}

/* --------------------------------------------------------------------------- */
/*    Static registration of the LAPACK routines                               */
/* --------------------------------------------------------------------------- */

bool OpenMPSCALAPACK::registerInterface()
{
    LAMAInterface& interface = LAMAInterfaceRegistry::getRegistry().modifyInterface( Context::Host );
    setInterface( interface.BLAS );
    return true;
}

/* --------------------------------------------------------------------------- */
/*    Static initialiazion at program start                                    */
/* --------------------------------------------------------------------------- */

bool OpenMPSCALAPACK::initialized = registerInterface();

} /* namespace lama */
