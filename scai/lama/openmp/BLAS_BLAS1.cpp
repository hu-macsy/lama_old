/**
 * @file BLAS_BLAS1.cpp
 *
 * @license
 * Copyright (c) 2009-2015
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
 * @brief Implementation of LAMA BLAS1 routines using wrappers to BLAS library
 * @author Lauretta Schubert
 * @date 05.07.2012
 * @since 1.0.0
 */

// hpp
#include <scai/lama/openmp/BLAS_BLAS1.hpp>

// local library
#include <scai/lama/openmp/BLASHelper.hpp>
#include <scai/lama/BLASInterface.hpp>
#include <scai/lama/LAMAInterfaceRegistry.hpp>
#include <scai/lama/cblas.hpp>

//Intel MPI need mpi.h to be included before stdio.h so this header comes first
#include <scai/tracing.hpp>

#include <scai/common/Settings.hpp>
#include <scai/common/macros/unused.hpp>

// boost
#include <boost/preprocessor.hpp>

// std
#include <cmath>

namespace scai
{

namespace lama
{

using tasking::SyncToken;
using common::getScalarType;

SCAI_LOG_DEF_LOGGER( BLAS_BLAS1::logger, "BLAS.BLAS1" )

/* ---------------------------------------------------------------------------------------*/
/*    scal wrapper for BLAS                                                               */
/* ---------------------------------------------------------------------------------------*/

template<typename ValueType>
static inline void wrapperScal( const int n, const ValueType alpha, ValueType* x, const int incX );

template<>
void wrapperScal( const int n, const float alpha, float* x, const int incX )
{
    cblas_sscal( n, alpha, x, incX );
}

template<>
void wrapperScal( const int n, const double alpha, double* x, const int incX )
{
    cblas_dscal( n, alpha, x, incX );
}

void wrapperScal( const int n, const ComplexFloat alpha, ComplexFloat* x, const int incX )
{
    // Attention: alpha is here passed by a pointer
    cblas_cscal( n, &alpha, x, incX );
}

template<>
void wrapperScal( const int n, const ComplexDouble alpha, ComplexDouble* x, const int incX )
{
    // Attention: alpha is here passed by a pointer
    cblas_zscal( n, &alpha, x, incX );
}

template<>
void wrapperScal( const int, const LongDouble, LongDouble*, const int )
{
    COMMON_THROWEXCEPTION( "LongDouble not supported by BLAS, please set USE_BLAS=0" )
}

template<>
void wrapperScal( const int, const ComplexLongDouble, ComplexLongDouble*, const int )
{
    COMMON_THROWEXCEPTION( "ComplexLongDouble not supported by BLAS, please set USE_BLAS=0" )
}

/* ---------------------------------------------------------------------------------------*/
/*    scal                                                                                */
/* ---------------------------------------------------------------------------------------*/

template<typename ValueType>
void BLAS_BLAS1::scal(
    const IndexType n,
    const ValueType alpha,
    ValueType* x,
    const IndexType incX,
    SyncToken* syncToken )
{
    if( incX <= 0 )
    {
        return;
    }

    SCAI_REGION( "BLAS.BLAS1.scal" )

    SCAI_LOG_DEBUG( logger,
                    "scal<" << getScalarType<ValueType>() << " n = " << n << ", alpha = " << alpha << ", x = " << x << ", incX = " << incX )

    if( syncToken )
    {
        SCAI_LOG_WARN( logger, "no asynchronous execution for openmp possible at this level." )
    }

    wrapperScal( n, alpha, x, incX );
}

/* ---------------------------------------------------------------------------------------*/
/*    nrm2                                                                                */
/* ---------------------------------------------------------------------------------------*/

template<typename ValueType>
static inline ValueType wrapperNrm2( const int n, const ValueType* x, const int incX );

template<>
float wrapperNrm2( const int n, const float* x, const int incX )
{
    return cblas_snrm2( n, x, incX );
}

template<>
double wrapperNrm2( const int n, const double* x, const int incX )
{
    return cblas_dnrm2( n, x, incX );
}

template<>
ComplexFloat wrapperNrm2( const int n, const ComplexFloat* x, const int incX )
{
    float res = cblas_scnrm2( n, x, incX );
    return ComplexFloat( res );
}

template<>
ComplexDouble wrapperNrm2( const int n, const ComplexDouble* x, const int incX )
{
    double res = cblas_dznrm2( n, x, incX );
    return ComplexDouble( res );
}

template<>
LongDouble wrapperNrm2( const int, const LongDouble*, const int )
{
    COMMON_THROWEXCEPTION( "LongDouble not supported by BLAS, please set USE_BLAS=0" )
    return static_cast<LongDouble>(0.0);
}

template<>
ComplexLongDouble wrapperNrm2( const int, const ComplexLongDouble*, const int )
{
    COMMON_THROWEXCEPTION( "LongDouble not supported by BLAS, please set USE_BLAS=0" )
    return static_cast<ComplexLongDouble>(0.0);
}

template<typename ValueType>
ValueType BLAS_BLAS1::nrm2( const IndexType n, const ValueType* x, const IndexType incX, SyncToken* syncToken )
{
    if( incX <= 0 )
    {
        return static_cast<ValueType>(0.0);
    }

    SCAI_REGION( "BLAS.BLAS1.nrm2" )

    SCAI_LOG_DEBUG( logger,
                    "nrm2<" << getScalarType<ValueType>() << ">, n = " << n << ", x = " << x << ", incX = " << incX )

    if( syncToken )
    {
        SCAI_LOG_WARN( logger, "no asynchronous execution for openmp possible at this level." )
    }

    return wrapperNrm2( n, x, incX );
}

/* ---------------------------------------------------------------------------------------*/
/*    asum                                                                                */
/* ---------------------------------------------------------------------------------------*/

template<typename ValueType>
static inline ValueType wrapperAsum( const int n, const ValueType* x, const int incX );

template<>
float wrapperAsum( const int n, const float* x, const int incX )
{
    return cblas_sasum( n, x, incX );
}

template<>
double wrapperAsum( const int n, const double* x, const int incX )
{
    return cblas_dasum( n, x, incX );
}

template<>
ComplexFloat wrapperAsum( const int n, const ComplexFloat* x, const int incX )
{
    float res = cblas_scasum( n, x, incX );
    return ComplexFloat( res );
}

template<>
ComplexDouble wrapperAsum( const int n, const ComplexDouble* x, const int incX )
{
    double res = cblas_dzasum( n, x, incX );
    return ComplexDouble( res );
}

template<>
LongDouble wrapperAsum( const int, const LongDouble*, const int )
{
    COMMON_THROWEXCEPTION( "LongDouble not supported by BLAS, please set USE_BLAS=0" )
    return static_cast<LongDouble>(0.0);
}

template<>
ComplexLongDouble wrapperAsum( const int, const ComplexLongDouble*, const int )
{
    COMMON_THROWEXCEPTION( "LongDouble not supported by BLAS, please set USE_BLAS=0" )
    return static_cast<ComplexLongDouble>(0.0);
}

template<typename ValueType>
ValueType BLAS_BLAS1::asum( const IndexType n, const ValueType* x, const IndexType incX, SyncToken* syncToken )
{
    SCAI_REGION( "BLAS.BLAS1.asum" )

    SCAI_LOG_DEBUG( logger,
                    "asum<" << getScalarType<ValueType>() << ">,  n = " << n << ", x = " << x << ", incX = " << incX )

    if( syncToken )
    {
        SCAI_LOG_WARN( logger, "asynchronous execution of BLAS1:asum not supported, do it synchronously." )
    }

    if( incX <= 0 )
    {
        return static_cast<ValueType>(0.0);
    }

    return wrapperAsum( n, x, incX );
}

/* ---------------------------------------------------------------------------------------*/
/*    iamax                                                                               */
/* ---------------------------------------------------------------------------------------*/

// Note: implicit conversion from int to IndexType and vice versa is done here
template<typename ValueType>
static IndexType wrapperIamax( const int n, const ValueType* x, const int incX );

template<>
IndexType wrapperIamax( const int n, const float* x, const int incX )
{
    return cblas_isamax( n, x, incX );
}

template<>
IndexType wrapperIamax( const int n, const double* x, const int incX )
{
    return cblas_idamax( n, x, incX );
}

template<>
IndexType wrapperIamax( const int n, const ComplexFloat* x, const int incX )
{
    return cblas_icamax( n, x, incX );
}

template<>
IndexType wrapperIamax( const int n, const ComplexDouble* x, const int incX )
{
    return cblas_izamax( n, x, incX );
}

template<>
IndexType wrapperIamax( const int, const LongDouble*, const int )
{
    COMMON_THROWEXCEPTION( "LongDouble not supported by BLAS, please set USE_BLAS=0" )
    return 0;
}

template<>
IndexType wrapperIamax( const int, const ComplexLongDouble*, const int )
{
    COMMON_THROWEXCEPTION( "LongDouble not supported by BLAS, please set USE_BLAS=0" )
    return 0;
}

template<typename ValueType>
IndexType BLAS_BLAS1::iamax( const IndexType n, const ValueType* x, const IndexType incX, SyncToken* syncToken )
{
    SCAI_REGION( "BLAS.BLAS1.iamax" )

    SCAI_LOG_INFO( logger,
                   "iamax<" << getScalarType<ValueType>() << ">, " << "n = " << n << ", x = " << x << ", incX = " << incX )

    if( syncToken )
    {
        SCAI_LOG_WARN( logger, "no asynchronous execution for openmp possible at this level." )
    }

    return wrapperIamax( n, x, incX );
}

/* ---------------------------------------------------------------------------------------*/
/*    swap                                                                                */
/* ---------------------------------------------------------------------------------------*/

template<typename ValueType>
static inline void wrapperSwap( const int n, ValueType* x, const int incX, ValueType* y, const int incY );

template<>
void wrapperSwap( const int n, float* x, const int incX, float* y, const int incY )
{
    cblas_sswap( n, x, incX, y, incY );
}

template<>
void wrapperSwap( const int n, double* x, const int incX, double* y, const int incY )
{
    cblas_dswap( n, x, incX, y, incY );
}

template<>
void wrapperSwap( const int n, ComplexFloat* x, const int incX, ComplexFloat* y, const int incY )
{
    cblas_cswap( n, x, incX, y, incY );
}

template<>
void wrapperSwap( const int n, ComplexDouble* x, const int incX, ComplexDouble* y, const int incY )
{
    cblas_zswap( n, x, incX, y, incY );
}

template<>
void wrapperSwap( const int, LongDouble*, const int, LongDouble*, const int )
{
    COMMON_THROWEXCEPTION( "LongDouble not supported by BLAS, please set USE_BLAS=0" )
}

template<>
void wrapperSwap( const int, ComplexLongDouble*, const int, ComplexLongDouble*, const int )
{
    COMMON_THROWEXCEPTION( "LongDouble not supported by BLAS, please set USE_BLAS=0" )
}

template<typename ValueType>
void BLAS_BLAS1::swap(
    const IndexType n,
    ValueType* x,
    const IndexType incX,
    ValueType* y,
    const IndexType incY,
    SyncToken* syncToken )
{
    SCAI_REGION( "BLAS.BLAS1.swap" )

    SCAI_LOG_DEBUG( logger,
                    "iamax<" << getScalarType<ValueType>() << ">" << ", n = " << n << ", x = " << x << ", incX = " << incX << ", y = " << y << ", incY = " << incY )

    if( ( incX <= 0 ) || ( incY <= 0 ) )
    {
        return;
    }

    if( syncToken )
    {
        SCAI_LOG_WARN( logger, "no asynchronous execution for openmp possible at this level." )
    }

    wrapperSwap( n, x, incX, y, incY );
}

/* ---------------------------------------------------------------------------------------*/
/*    copy                                                                                */
/* ---------------------------------------------------------------------------------------*/

template<typename ValueType>
static inline void wrapperCopy( const int n, const ValueType* x, const int incX, ValueType* y, const int incY );

template<>
void wrapperCopy( const int n, const float* x, const int incX, float* y, const int incY )
{
    cblas_scopy( n, x, incX, y, incY );
}

template<>
void wrapperCopy( const int n, const double* x, const int incX, double* y, const int incY )
{
    cblas_dcopy( n, x, incX, y, incY );
}

template<>
void wrapperCopy( const int n, const ComplexFloat* x, const int incX, ComplexFloat* y, const int incY )
{
    cblas_ccopy( n, x, incX, y, incY );
}

template<>
void wrapperCopy( const int n, const ComplexDouble* x, const int incX, ComplexDouble* y, const int incY )
{
    cblas_zcopy( n, x, incX, y, incY );
}

template<>
void wrapperCopy( const int, const LongDouble*, const int, LongDouble*, const int )
{
    COMMON_THROWEXCEPTION( "LongDouble not supported by BLAS, please set USE_BLAS=0" )
}


template<>
void wrapperCopy( const int, const ComplexLongDouble*, const int, ComplexLongDouble*, const int )
{
    COMMON_THROWEXCEPTION( "LongDouble not supported by BLAS, please set USE_BLAS=0" )
}

template<typename ValueType>
void BLAS_BLAS1::copy(
    const IndexType n,
    const ValueType* x,
    const IndexType incX,
    ValueType* y,
    const IndexType incY,
    SyncToken* syncToken )
{
    SCAI_REGION( "BLAS.BLAS1.copy" )

    SCAI_LOG_DEBUG( logger,
                    "copy<" << getScalarType<ValueType>() << ">, " << "n = " << n << ", x = " << x << ", incX = " << incX << ", y = " << y << ", incY = " << incY )

    if( ( incX <= 0 ) || ( incY <= 0 ) )
    {
        return;
    }

    if( syncToken )
    {
        SCAI_LOG_WARN( logger, "no asynchronous execution for openmp possible at this level." )
    }

    // wrapper: uses overloading and converts arguments if required

    wrapperCopy( n, x, incX, y, incY );
}

/* ---------------------------------------------------------------------------------------*/
/*    axpy                                                                                */
/* ---------------------------------------------------------------------------------------*/

template<typename ValueType>
static inline void wrapperAxpy(
    const int n,
    ValueType alpha,
    const ValueType* x,
    const int incX,
    ValueType* y,
    const int incY );

template<>
void wrapperAxpy( const int n, float alpha, const float* x, const int incX, float* y, const int incY )
{
    cblas_saxpy( n, alpha, x, incX, y, incY );
}

template<>
void wrapperAxpy( const int n, double alpha, const double* x, const int incX, double* y, const int incY )
{
    cblas_daxpy( n, alpha, x, incX, y, incY );
}

template<>
void wrapperAxpy(
    const int n,
    ComplexFloat alpha,
    const ComplexFloat* x,
    const int incX,
    ComplexFloat* y,
    const int incY )
{
    // Attention: alpha is here passed by a pointer
    cblas_caxpy( n, &alpha, x, incX, y, incY );
}

template<>
void wrapperAxpy(
    const int n,
    ComplexDouble alpha,
    const ComplexDouble* x,
    const int incX,
    ComplexDouble* y,
    const int incY )
{
    // Attention: alpha is here passed by a pointer
    cblas_zaxpy( n, &alpha, x, incX, y, incY );
}

template<>
void wrapperAxpy( const int, LongDouble, const LongDouble*, const int, LongDouble*, const int )
{
    COMMON_THROWEXCEPTION( "LongDouble not supported by BLAS, please set USE_BLAS=0" )
}

template<>
void wrapperAxpy( const int, ComplexLongDouble, const ComplexLongDouble*, const int, ComplexLongDouble*, const int )
{
    COMMON_THROWEXCEPTION( "LongDouble not supported by BLAS, please set USE_BLAS=0" )
}

template<typename ValueType>
void BLAS_BLAS1::axpy(
    const IndexType n,
    const ValueType alpha,
    const ValueType* x,
    const IndexType incX,
    ValueType* y,
    const IndexType incY,
    SyncToken* syncToken )
{
    SCAI_REGION( "BLAS.BLAS1.axpy" )

    SCAI_LOG_DEBUG( logger,
                    "axpy<" << getScalarType<ValueType>() << ">, " << "n = " << n << ", alpha = " << alpha << ", x = " << x << ", incX = " << incX << ", y = " << y << ", incY = " << incY )

    if( ( incX <= 0 ) || ( incY <= 0 ) )
    {
        return;
    }

    if( syncToken )
    {
        SCAI_LOG_WARN( logger, "no asynchronous execution for openmp possible at this level." )
    }

    wrapperAxpy( n, alpha, x, incX, y, incY );
}

/* ---------------------------------------------------------------------------------------*/
/*    dot                                                                                 */
/* ---------------------------------------------------------------------------------------*/

template<typename ValueType>
static inline ValueType wrapperDot(
    const int n,
    const ValueType* x,
    const int incX,
    const ValueType* y,
    const int incY );

template<>
float wrapperDot( const int n, const float* x, const int incX, const float* y, const int incY )
{
    return cblas_sdot( n, x, incX, y, incY );
}

template<>
double wrapperDot( const int n, const double* x, const int incX, const double* y, const int incY )
{
    return cblas_ddot( n, x, incX, y, incY );
}

template<>
ComplexFloat wrapperDot( const int n, const ComplexFloat* x, const int incX, const ComplexFloat* y, const int incY )
{
    ComplexFloat dotu;
    cblas_cdotu_sub( n, x, incX, y, incY, &dotu );
    return dotu;
}

template<>
ComplexDouble wrapperDot( const int n, const ComplexDouble* x, const int incX, const ComplexDouble* y, const int incY )
{
    ComplexDouble dotu;
    cblas_zdotu_sub( n, x, incX, y, incY, &dotu );
    return dotu;
}

template<>
LongDouble wrapperDot( const int, const LongDouble*, const int, const LongDouble*, const int )
{
    COMMON_THROWEXCEPTION( "LongDouble not supported by BLAS, please set USE_BLAS=0" )
    return static_cast<LongDouble>(0.0);
}

template<>
ComplexLongDouble wrapperDot( const int, const ComplexLongDouble*, const int, const ComplexLongDouble*, const int )
{
    COMMON_THROWEXCEPTION( "LongDouble not supported by BLAS, please set USE_BLAS=0" )
    return static_cast<ComplexLongDouble>(0.0);
}

template<typename ValueType>
ValueType BLAS_BLAS1::dot(
    const IndexType n,
    const ValueType* x,
    const IndexType incX,
    const ValueType* y,
    const IndexType incY,
    SyncToken* syncToken )
{
    SCAI_REGION( "BLAS.BLAS1.dot" )

    SCAI_LOG_DEBUG( logger,
                    "dot<" << getScalarType<ValueType>() << ">, " << "n = " << n << ", x = " << x << ", incX = " << incX << ", y = " << y << ", incY = " << incY )

    if( ( incX <= 0 ) || ( incY <= 0 ) )
    {
        return static_cast<ValueType>(0.0);
    }

    if( syncToken )
    {
        SCAI_LOG_WARN( logger, "no asynchronous execution for openmp possible at this level." )
    }

    return wrapperDot( static_cast<int>( n ), x, static_cast<int>( incX ), y, static_cast<int>( incY ) );
}

/* --------------------------------------------------------------------------- */
/*     Template instantiations via registration routine                        */
/* --------------------------------------------------------------------------- */

void BLAS_BLAS1::setInterface( BLASInterface& BLAS )
{
    bool useBLAS = false;
    int level = 0;

    useBLAS = common::Settings::getEnvironment( level, "USE_BLAS" );

    if( !useBLAS || ( level <= 0 ) )
    {
        SCAI_LOG_INFO( logger, "BLAS1 wrapper routines for Host Interface are disabled (USE_BLAS not set or 0)" )
        return;
    }
    else if( level > 1 )
    {
        // only level 2 or level 3 wrappers might be used
        SCAI_LOG_INFO( logger,
                       "BLAS1 wrapper routines for Host Interface are disabled (USE_BLAS = " << level << ")" )
        return;
    }

    SCAI_LOG_INFO( logger, "set BLAS1 wrapper routines for Host Context in Interface" )

    // REGISTER1: give these routines priority in case of overriding

    // Note: macro takes advantage of same name for routines and type definitions
    //       ( e.g. routine CUDABLAS1::sum<ValueType> is set for BLAS::BLAS1::sum variable

#define LAMA_BLAS1_REGISTER(z, I, _)                                          \
    LAMA_INTERFACE_REGISTER1_T( BLAS, scal, ARITHMETIC_HOST_TYPE_##I )        \
    LAMA_INTERFACE_REGISTER1_T( BLAS, nrm2, ARITHMETIC_HOST_TYPE_##I )        \
    LAMA_INTERFACE_REGISTER1_T( BLAS, asum, ARITHMETIC_HOST_TYPE_##I )        \
    LAMA_INTERFACE_REGISTER1_T( BLAS, iamax, ARITHMETIC_HOST_TYPE_##I )       \
    LAMA_INTERFACE_REGISTER1_T( BLAS, swap, ARITHMETIC_HOST_TYPE_##I )        \
    LAMA_INTERFACE_REGISTER1_T( BLAS, copy, ARITHMETIC_HOST_TYPE_##I )        \
    LAMA_INTERFACE_REGISTER1_T( BLAS, axpy, ARITHMETIC_HOST_TYPE_##I )        \
    LAMA_INTERFACE_REGISTER1_T( BLAS, dot, ARITHMETIC_HOST_TYPE_##I )         \

    BOOST_PP_REPEAT( ARITHMETIC_HOST_EXT_TYPE_CNT, LAMA_BLAS1_REGISTER, _ )

#undef LAMA_BLAS1_REGISTER
}

/* --------------------------------------------------------------------------- */
/*    Static registration of the BLAS1 routines                                */
/* --------------------------------------------------------------------------- */

bool BLAS_BLAS1::registerInterface()
{
    LAMAInterface& interface = LAMAInterfaceRegistry::getRegistry().modifyInterface( hmemo::context::Host );
    setInterface( interface.BLAS );
    return true;
}

/* --------------------------------------------------------------------------- */
/*    Static initialiazion at program start                                    */
/* --------------------------------------------------------------------------- */

bool BLAS_BLAS1::initialized = registerInterface();

} /* end namespace lama */

} /* end namespace scai */
