/**
 * @file SCALAPACKHelper.hpp
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
 * @brief SCALAPACKHelper.hpp
 * @author lschubert
 * @date 02.07.2012
 * $Id$
 */
#ifndef LAMA_SCALAPACKHELPER_HPP_
#define LAMA_SCALAPACKHELPER_HPP_

// macros
#include <lama/macros/unused.hpp>

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

#ifdef __cplusplus
extern "C"
{
#endif

    extern void F77_descinit( int*, int*, int*, int*, int*, int*, int*, int*, int*, int* );
    extern int F77_numroc( int*, int*, int*, int*, int* );

#ifdef __cplusplus
}
#endif

struct SCALAPACKHelper
{
public:

    static void BLACS_Context( int* context );

    static void GridInit( int* context, const char* kind, const int npRow, const int npCol );

    static void GridInfo( const int Context, const int npRow, const int npCol, int* myRow, int* myCol );

    static void GridExit( const int context );

    static void BlacsFreeBuf( const int context );

    static int DescInit( int* descA, int m, int n, int mb, int nb, int irscr, int icsrc, int ictxt, int lld );

    static int numRoC( int n, int nb, int iproc, int nprocs );

    /**
     *  ILCM computes and returns the Least Common Multiple (LCM) of two
     *  positive integers M and N. In fact the routine computes the greatest
     *  common divisor (GCD) and use the fact that M*N = GCD*LCM.
     */
    static int ILCM( int m, int n );

};

#ifdef LAMA_MKL_SCALAPACK
void SCALAPACKHelper::BLACS_Context( int* context )
{
    MKL_INT what = 0; // ask for system context
    MKL_INT val = -1;
    F77_blacs_get( &val, &what, context );
}
#else
void SCALAPACKHelper::BLACS_Context( int* UNUSED(context) )
{
    // ERROR
    printf( "ERROR: BLACS not available\n" );
}
#endif

#ifdef LAMA_MKL_SCALAPACK
void SCALAPACKHelper::GridInit( int* context, const char* kind, const int npRow, const int npCol )
{
    F77_blacs_gridinit( context, (char*) kind, (int*) &npRow, (int*) &npCol );
}
#else
void SCALAPACKHelper::GridInit(
    int* UNUSED(context),
    const char* UNUSED(kind),
    const int UNUSED(npRow),
    const int UNUSED(npCol) )
{
}
#endif

#ifdef LAMA_MKL_SCALAPACK
void SCALAPACKHelper::GridInfo( const int Context, const int npRow, const int npCol, int* myRow, int* myCol )
{
    F77_blacs_gridinfo( (int*) &Context, (int*) &npRow, (int*) &npCol, myRow, myCol );
}
#else
void SCALAPACKHelper::GridInfo(
    const int UNUSED(Context),
    const int UNUSED(npRow),
    const int UNUSED(npCol),
    int* UNUSED(myRow),
    int* UNUSED(myCol) )
{
}
#endif

#ifdef LAMA_MKL_SCALAPACK
void SCALAPACKHelper::GridExit( const int context )
{
    F77_blacs_gridexit( (int*)( &context ) );
}
#else
void SCALAPACKHelper::GridExit( const int UNUSED(context) )
{
}
#endif

#ifdef LAMA_MKL_SCALAPACK
void SCALAPACKHelper::BlacsFreeBuf( const int context )
{
    MKL_INT i_one = 1;
    F77_blacs_freebuff( (int*)( &context ), &i_one );
}
#else
void SCALAPACKHelper::BlacsFreeBuf( const int UNUSED(context) )
{
}
#endif

#ifdef LAMA_MKL_SCALAPACK
int SCALAPACKHelper::DescInit( int* descA, int m, int n, int mb, int nb, int irscr, int icsrc, int ictxt, int lld)
{
    int info = 0;
    F77_descinit( descA, &m, &n, &mb, &nb, &irscr, &icsrc, &ictxt, &lld, &info );
    return info;
}
#else
int SCALAPACKHelper::DescInit(
    int* UNUSED(descA),
    int UNUSED(m),
    int UNUSED(n),
    int UNUSED(mb),
    int UNUSED(nb),
    int UNUSED(irscr),
    int UNUSED(icsrc),
    int UNUSED(ictxt),
    int UNUSED(lld) )
{
    return 0;
}
#endif

#ifdef LAMA_MKL_SCALAPACK
int SCALAPACKHelper::numRoC( int n, int nb, int iproc, int nprocs )
{
    int iSrcProc = 0; // first row/col always on first processor
    return F77_numroc (&n, &nb, &iproc, &iSrcProc, &nprocs );
}
#else
int SCALAPACKHelper::numRoC( int UNUSED(n), int UNUSED(nb), int UNUSED(iproc), int UNUSED(nprocs) )
{
    return -1;
}
#endif

int SCALAPACKHelper::ILCM( int m, int n )
{
    int ia, iq, ir;
    int ilcm;

    if( m >= n )
    {
        ia = m;
        ilcm = n;
    }
    else
    {
        ia = n;
        ilcm = m;
    }

    iq = ia / ilcm;
    ir = ia - iq * ilcm;
    while( ir != 0 )
    {
        ia = ilcm;
        ilcm = ir;
        iq = ia / ilcm;
        ir = ia - iq * ilcm;
    }
    return ( m * n ) / ilcm;
}

#endif // LAMA_SCALAPACKHELPER_HPP_
