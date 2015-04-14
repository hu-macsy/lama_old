/**
 * @file CUDATexVector.hpp
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
 * @brief Routines for indexed read access of a vector X
 * @author Thomas Brandes, Jiri Kraus
 * @date 04.07.2012
 * @since 1.0.0
 */

/* --------------------------------------------------------------------------- */

/************************************************************************************
 *
 *  Global variable for vector texture is needed for all supported arithmetic types
 *
 *  Due to static declaration it is safe to use it in different source files
 *
 *  Some algorihtms use two vectors at the same type (might be unsupported on some devices)
 *
 *************************************************************************************/

static texture<int4,1> texVectorZXref;

static texture<float2,1> texVectorCXref;

static texture<float,1> texVectorSXref;

static texture<int2,1> texVectorDXref;

static texture<int,1> texVectorIref;

__inline__ static void vectorBindTexture( const float* vector )
{
    LAMA_CUDA_RT_CALL( cudaBindTexture( NULL, texVectorSXref, vector ), "bind float vector x to texture" )
}

__inline__ static void vectorBindTexture( const double* vector )
{
    LAMA_CUDA_RT_CALL( cudaBindTexture( NULL, texVectorDXref, vector ), "bind double vector x to texture" )
}

__inline__ static void vectorBindTexture( const ComplexFloat* vector )
{
    LAMA_CUDA_RT_CALL( cudaBindTexture( NULL, texVectorCXref, vector ), "bind ComplexFloat vector x to texture" )
}

__inline__ static void vectorBindTexture( const ComplexDouble* vector )
{
    LAMA_CUDA_RT_CALL( cudaBindTexture( NULL, texVectorZXref, vector ), "bind ComplexDouble vector x to texture" )
}

__inline__ static void vectorBindTexture( const int* vector )
{
    LAMA_CUDA_RT_CALL( cudaBindTexture( NULL, texVectorIref, vector ), "bind int vector x to texture" )
}

__inline__ static void vectorUnbindTexture( const float* )
{
    LAMA_CUDA_RT_CALL( cudaUnbindTexture( texVectorSXref ), "unbind float vector x from texture" )
}

__inline__ static void vectorUnbindTexture( const double* )
{
    LAMA_CUDA_RT_CALL( cudaUnbindTexture( texVectorDXref ), "unbind double vector x from texture" )
}

__inline__ static void vectorUnbindTexture( const ComplexFloat* )
{
    LAMA_CUDA_RT_CALL( cudaUnbindTexture( texVectorCXref ), "unbind ComplexFloat vector x from texture" )
}

__inline__ static void vectorUnbindTexture( const ComplexDouble* )
{
    LAMA_CUDA_RT_CALL( cudaUnbindTexture( texVectorZXref ), "unbind ComplexDouble vector x from texture" )
}

__inline__ static void vectorUnbindTexture( const int* )
{
    LAMA_CUDA_RT_CALL( cudaUnbindTexture( texVectorIref ), "unbind int vector x from texture" )
}

template<typename ValueType, bool useTexture>
__inline__ __device__
static ValueType fetchVectorX( const ValueType* const x, const int i )
{
    return x[i];
}

template<>
__inline__ __device__
float fetchVectorX<float, true>( const float* const, const int i )
{
    return tex1Dfetch( texVectorSXref, i );
}

template<>
__inline__ __device__
double fetchVectorX<double, true>( const double* const, const int i )
{
    int2 v = tex1Dfetch( texVectorDXref, i );
    return __hiloint2double( v.y, v.x );
}

template<>
__inline__ __device__
int fetchVectorX<int, true>( const int* const, const int i )
{
    return tex1Dfetch( texVectorIref, i );
}

template<>
__inline__ __device__
ComplexFloat fetchVectorX<ComplexFloat, true>( const ComplexFloat* const, const int i )
{
    float2 v = tex1Dfetch( texVectorCXref, i );
    return ComplexFloat(v.x, v.y);
}

template<>
__inline__ __device__
ComplexDouble fetchVectorX<ComplexDouble, true>( const ComplexDouble* const, const int i )
{
    int4 u = tex1Dfetch( texVectorZXref, i );
    return ComplexDouble( __hiloint2double( u.y, u.x ), __hiloint2double( u.w, u.z));
}

