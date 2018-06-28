/**
 * @file CUDATexVector.hpp
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
 * @brief Routines for indexed read access of a vector X
 * @author Thomas Brandes, Jiri Kraus
 * @date 04.07.2012
 */

#pragma once

#include <scai/common/cuda/CUDAError.hpp>

namespace scai
{

/************************************************************************************
 *
 *  Global variable for vector texture is needed for all supported arithmetic types
 *
 *  Due to static declaration it is safe to use it in different source files
 *
 *  Some algorihtms use two vectors at the same type (might be unsupported on some devices)
 *
 *************************************************************************************/

static texture<int4, 1> texVectorZXref;

static texture<float2, 1> texVectorCXref;

static texture<float, 1> texVectorSXref;

static texture<int2, 1> texVectorDXref;

static texture<int, 1> texVectorIref;

static texture<long, 1> texVectorLref;

static texture<unsigned int, 1> texVectorUref;

static texture<unsigned long, 1> texVectorULref;

/* --------------------------------------------------------------------------- */
/*  Bind and Unbind device array data to the texture                           */
/* --------------------------------------------------------------------------- */

__inline__ static void vectorBindTexture( const float* vector )
{
    SCAI_CUDA_RT_CALL( cudaBindTexture( NULL, texVectorSXref, vector ), "bind float vector x to texture" )
}

__inline__ static void vectorBindTexture( const double* vector )
{
    SCAI_CUDA_RT_CALL( cudaBindTexture( NULL, texVectorDXref, vector ), "bind double vector x to texture" )
}

#ifdef SCAI_COMPLEX_SUPPORTED

__inline__ static void vectorBindTexture( const ComplexFloat* vector )
{
    SCAI_CUDA_RT_CALL( cudaBindTexture( NULL, texVectorCXref, vector ), "bind ComplexFloat vector x to texture" )
}

__inline__ static void vectorBindTexture( const ComplexDouble* vector )
{
    SCAI_CUDA_RT_CALL( cudaBindTexture( NULL, texVectorZXref, vector ), "bind ComplexDouble vector x to texture" )
}

#endif

__inline__ static void vectorBindTexture( const int* vector )
{
    SCAI_CUDA_RT_CALL( cudaBindTexture( NULL, texVectorIref, vector ), "bind int vector x to texture" )
}

__inline__ static void vectorBindTexture( const long* vector )
{
    SCAI_CUDA_RT_CALL( cudaBindTexture( NULL, texVectorLref, vector ), "bind int vector x to texture" )
}

__inline__ static void vectorBindTexture( const unsigned int* vector )
{
    SCAI_CUDA_RT_CALL( cudaBindTexture( NULL, texVectorUref, vector ), "bind unsigned int vector x to texture" )
}

__inline__ static void vectorBindTexture( const unsigned long* vector )
{
    SCAI_CUDA_RT_CALL( cudaBindTexture( NULL, texVectorULref, vector ), "bind unsigned long vector x to texture" )
}

__inline__ static void vectorUnbindTexture( const float* )
{
    SCAI_CUDA_RT_CALL( cudaUnbindTexture( texVectorSXref ), "unbind float vector x from texture" )
}

__inline__ static void vectorUnbindTexture( const double* )
{
    SCAI_CUDA_RT_CALL( cudaUnbindTexture( texVectorDXref ), "unbind double vector x from texture" )
}

#ifdef SCAI_COMPLEX_SUPPORTED

__inline__ static void vectorUnbindTexture( const ComplexFloat* )
{
    SCAI_CUDA_RT_CALL( cudaUnbindTexture( texVectorCXref ), "unbind ComplexFloat vector x from texture" )
}

__inline__ static void vectorUnbindTexture( const ComplexDouble* )
{
    SCAI_CUDA_RT_CALL( cudaUnbindTexture( texVectorZXref ), "unbind ComplexDouble vector x from texture" )
}

#endif

__inline__ static void vectorUnbindTexture( const int* )
{
    SCAI_CUDA_RT_CALL( cudaUnbindTexture( texVectorIref ), "unbind int vector x from texture" )
}

__inline__ static void vectorUnbindTexture( const unsigned int* )
{
    SCAI_CUDA_RT_CALL( cudaUnbindTexture( texVectorUref ), "unbind unsigned int vector x from texture" )
}

__inline__ static void vectorUnbindTexture( const long* )
{
    SCAI_CUDA_RT_CALL( cudaUnbindTexture( texVectorLref ), "unbind long vector x from texture" )
}

__inline__ static void vectorUnbindTexture( const unsigned long* )
{
    SCAI_CUDA_RT_CALL( cudaUnbindTexture( texVectorULref ), "unbind unsigned long vector x from texture" )
}

/* --------------------------------------------------------------------------- */
/*  GPU device routines to access values of the texture                        */
/* --------------------------------------------------------------------------- */

template<typename ValueType, bool useTexture>
__inline__ __device__
static ValueType fetchVectorX( const ValueType* const x, const IndexType i )
{
    return x[i];
}

template<>
__inline__ __device__
float fetchVectorX<float, true>( const float* const, const IndexType i )
{
    return tex1Dfetch( texVectorSXref, i );
}

template<>
__inline__ __device__
double fetchVectorX<double, true>( const double* const, const IndexType i )
{
    int2 v = tex1Dfetch( texVectorDXref, i );
    return __hiloint2double( v.y, v.x );
}

template<>
__inline__ __device__
int fetchVectorX<int, true>( const int* const, const IndexType i )
{
    return tex1Dfetch( texVectorIref, i );
}

template<>
__inline__ __device__
unsigned int fetchVectorX<unsigned int, true>( const unsigned int* const, const IndexType i )
{
    return tex1Dfetch( texVectorUref, i );
}

#ifdef SCAI_COMPLEX_SUPPORTED

template<>
__inline__ __device__
ComplexFloat fetchVectorX<ComplexFloat, true>( const ComplexFloat* const, const IndexType i )
{
    float2 v = tex1Dfetch( texVectorCXref, i );
    return ComplexFloat( v.x, v.y );
}

template<>
__inline__ __device__
ComplexDouble fetchVectorX<ComplexDouble, true>( const ComplexDouble* const, const IndexType i )
{
    int4 u = tex1Dfetch( texVectorZXref, i );
    return ComplexDouble( __hiloint2double( u.y, u.x ), __hiloint2double( u.w, u.z ) );
}

#endif

}
