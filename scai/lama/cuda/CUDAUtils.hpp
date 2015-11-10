/**
 * @file CUDAUtils.hpp
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
 * @brief Implementation of general utilities with CUDA
 * @author Thomas Brandes
 * @date 02.07.2012
 * @since 1.0.0
 */

#pragma once

// for dll_import
#include <scai/common/config.hpp>

// internal scai libraries
#include <scai/logging.hpp>

#include <scai/common/SCAITypes.hpp>
#include <scai/common/Assert.hpp>

namespace scai
{

namespace lama
{

/** General utilities of the LAMA Interface implemented in CUDA  */

class COMMON_DLL_IMPORTEXPORT CUDAUtils
{
public:

    /*  CUDA implementation of UtilKernelTrait::validIndexes  */

    static bool validIndexes( const IndexType array[], const IndexType n, const IndexType size );

    /*  CUDA implementation of UtilKernelTrait::sum  */

    template<typename ValueType>
    static ValueType sum( const ValueType array[], const IndexType n );

    /*  CUDA implementation of UtilKernelTrait::setVal  */

    template<typename ValueType>
    static void setVal( ValueType array[], const IndexType n, const ValueType val );

    /*  CUDA implementation of UtilKernelTrait::setOrder  */

    template<typename ValueType>
    static void setOrder( ValueType array[], const IndexType n );

    /*  CUDA implementation of UtilKernelTrait::getValue  */

    template<typename ValueType>
    static ValueType getValue( const ValueType* array, const IndexType i );

    /** CUDA implementation for UtilKernelTrait::scale. */

    template<typename ValueType>
    static void scale( ValueType values[], const ValueType value, const IndexType n );

    /** CUDA implementation for UtilKernelTrait::setScale. */

    template<typename ValueType,typename otherValueType>
    static void setScale(
        ValueType outValues[],
        const ValueType value,
        const otherValueType inValues[],
        const IndexType n );

    /** CUDA function implements UtilKernelTrait::maxval */

    template<typename ValueType>
    static ValueType maxval( const ValueType array[], const IndexType n );

    /** CUDA function implements UtilKernelTrait::absMaxVal */

    template<typename ValueType>
    static ValueType absMaxVal( const ValueType array[], const IndexType n );

    /** CUDA function implements UtilKernelTrait::absMaxDiffVal */

    template<typename ValueType>
    static ValueType absMaxDiffVal( const ValueType array1[], const ValueType array2[], const IndexType n );

    /** CUDA implementation for UtilKernelTrait::isSorted */

    template<typename ValueType>
    static bool isSorted( const ValueType array[], const IndexType n, bool acending );

    /** CUDA implementation for UtilKernelTrait::set */

    template<typename ValueType,typename otherValueType>
    static void set( ValueType out[], const otherValueType in[], const IndexType n );

    /** CUDA implementation for UtilKernelTrait::setGather, out[i]] = in[ indexes[i] ] */

    template<typename ValueType,typename otherValueType>
    static void setGather( ValueType out[], const otherValueType in[], const IndexType indexes[], const IndexType n );

    /** CUDA implementation for UtilKernelTrait::setScatter, out[ indexes[i] ] = in [i] */

    template<typename ValueType,typename otherValueType>
    static void setScatter( ValueType out[], const IndexType indexes[], const otherValueType in[], const IndexType n );

    /** CUDA implementation for UtilKernelTrait::invert */

    template<typename ValueType>
    static void invert( ValueType array[], const IndexType n );

private:

    SCAI_LOG_DECL_STATIC_LOGGER( logger )

    /** Routine that registers all methods at the kernel registry. */

    static void registerKernels( bool deleteFlag );

    /** Constructor for registration. */

    CUDAUtils();

    /** Destructor for unregistration. */

    ~CUDAUtils();

    /** Static variable for registration at static initialization. */

    static CUDAUtils guard;
};

} /* end namespace lama */

} /* end namespace scai */
