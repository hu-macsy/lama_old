/**
 * @file sparsekernel/mic/MICMKLCSRWrapper.hpp
 *
 * @license
 * Copyright (c) 2009-2016
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
 * @brief ToDo: Missing description in ./sparsekernel/mic/MICMKLCSRWrapper.hpp
 * @author Eric Schricker
 * @date 09.02.2016
 */

#pragma once

// internal scai libraries
#include <scai/sparsekernel/mic/MICMKLCSRTrait.hpp>

// extern
#include <mkl_spblas.h>

namespace scai
{

namespace sparsekernel
{

template<typename ValueType>
class COMMON_DLL_IMPORTEXPORT MICMKLCSRWrapper;

#define MICMKLCSRWRAPPER_DEF( ValueType, MKLCSRValueType, prefix ) 								                                \
template<>																												        \
class COMMON_DLL_IMPORTEXPORT MICMKLCSRWrapper<ValueType>								                              			\
{																														        \
public:																													        \
	  typedef MICMKLCSRTrait::BLASIndexType BLASIndexType;																	    \
	  typedef MICMKLCSRTrait::BLASTrans BLASTrans;																	            \
	  typedef MICMKLCSRTrait::BLASMatrix BLASMatrix;																	        \
																														        \
    static BLASIndexType csr2csc(                                                                                               \
        const BLASIndexType* job,                                                                                               \
        const BLASIndexType n,                                                                                                  \
        const ValueType *Acsr,                                                                                                  \
        const BLASIndexType *AJ0,                                                                                               \
        const BLASIndexType *AI0,                                                                                               \
        ValueType *Acsc,                                                                                                        \
        BLASIndexType *AJ1,                                                                                                     \
        BLASIndexType *AI1 )                                                                                                    \
    {                                                                                                                           \
	      BLASIndexType info;                                                                                                   \
        MICMKLCSR_BLAS_NAME( csrcsc, prefix )( const_cast<BLASIndexType*>( job ), const_cast<BLASIndexType*>( &n ),             \
                                            const_cast<MKLCSRValueType*>( reinterpret_cast<const MKLCSRValueType*>( Acsr ) ),   \
                                            const_cast<BLASIndexType*>( AJ0 ), const_cast<BLASIndexType*>( AI0 ),               \
                                            const_cast<MKLCSRValueType*>( reinterpret_cast<MKLCSRValueType*>( Acsc ) ),         \
                                            const_cast<BLASIndexType*>( AJ1 ), const_cast<BLASIndexType*>( AI1 ), &info );      \
        return info;                                                                                                            \
    }                                                                                                                           \
																														        \
    static void csrmv(                                                                                                          \
        const BLASTrans transA,                                                                                                 \
        const BLASIndexType m,                                                                                                  \
        const BLASIndexType k,                                                                                                  \
        const ValueType alpha,                                                                                                  \
        const BLASMatrix descrA,                                                                                                \
        const ValueType *val,                                                                                                   \
        const BLASIndexType *indx,                                                                                              \
        const BLASIndexType *pntrb,                                                                                             \
        const BLASIndexType *pntre,                                                                                             \
        const ValueType *x,                                                                                                     \
        const ValueType beta,                                                                                                   \
        ValueType *y )                                                                                                          \
    {                                                                                                                           \
        MICMKLCSR_BLAS_NAME( csrmv, prefix )( const_cast<BLASTrans*>( &transA ), const_cast<BLASIndexType*>( &m ),              \
                                           const_cast<BLASIndexType*>( &k ),                                                    \
                                           const_cast<MKLCSRValueType*>( reinterpret_cast<const MKLCSRValueType*>( &alpha ) ),  \
                                           const_cast<char*>( descrA ),                                                         \
                                           const_cast<MKLCSRValueType*>( reinterpret_cast<const MKLCSRValueType*>( val ) ),     \
                                           const_cast<BLASIndexType*>( indx ), const_cast<BLASIndexType*>( pntrb ),             \
                                           const_cast<BLASIndexType*>( pntre ),                                                 \
                                           const_cast<MKLCSRValueType*>( reinterpret_cast<const MKLCSRValueType*>( x ) ),       \
                                           const_cast<MKLCSRValueType*>( reinterpret_cast<const MKLCSRValueType*>( &beta ) ),   \
                                           reinterpret_cast<MKLCSRValueType*>( y ) );                                           \
    }                                                                                                                           \
};

MICMKLCSRWRAPPER_DEF( float, float, s )
MICMKLCSRWRAPPER_DEF( double, double, d )

#ifdef SCAI_COMPLEX_SUPPORTED
    MICMKLCSRWRAPPER_DEF( ComplexFloat,  MKL_Complex8,  c )
    MICMKLCSRWRAPPER_DEF( ComplexDouble, MKL_Complex16, z )
#endif

#undef MICMKLCSRWRAPPER_DEF

} /* end namespace sparsekernel */

} /* end namespace scai */
