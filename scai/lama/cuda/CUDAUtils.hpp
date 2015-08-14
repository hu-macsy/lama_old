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

// others
#include <scai/lama/LAMATypes.hpp>

#include <scai/lama/exception/LAMAAssert.hpp>

// logging
#include <scai/logging.hpp>

namespace scai
{

namespace lama
{

/** General utilities of the LAMA Interface implemented in CUDA  */

class COMMON_DLL_IMPORTEXPORT CUDAUtils
{
public:

    /*  CUDA implementation of UtilsInterface::validIndexes  */

    static bool validIndexes( const IndexType array[], const IndexType n, const IndexType size );

    template<typename ValueType>
    static ValueType sum( const ValueType array[], const IndexType n );

    template<typename ValueType>
    static void setVal( ValueType array[], const IndexType n, const ValueType val );

    template<typename ValueType>
    static void setOrder( ValueType array[], const IndexType n );

    template<typename ValueType>
    static ValueType getValue( const ValueType* array, const IndexType i );

    /** CUDA implementation for UtilsInterface::Transform::scale. */

    template<typename ValueType>
    static void scale( ValueType values[], const ValueType value, const IndexType n );

    /** CUDA implementation for UtilsInterface::Copy::setScale. */

    template<typename ValueType,typename otherValueType>
    static void setScale(
        ValueType outValues[],
        const ValueType value,
        const otherValueType inValues[],
        const IndexType n );

    /** CUDA function implements UtilsInterface::Reductions::maxval */

    template<typename ValueType>
    static ValueType maxval( const ValueType array[], const IndexType n );

    /** CUDA function implements UtilsInterface::Reductions::absMaxVal */

    template<typename ValueType>
    static ValueType absMaxVal( const ValueType array[], const IndexType n );

    /** CUDA function implements UtilsInterface::Reductions::absMaxDiffVal */

    template<typename ValueType>
    static ValueType absMaxDiffVal( const ValueType array1[], const ValueType array2[], const IndexType n );

    /** CUDA implementation for UtilsInterface::Reductions::isSorted */

    template<typename ValueType>
    static bool isSorted( const ValueType array[], const IndexType n, bool acending );

    template<typename ValueType,typename otherValueType>
    static void set( ValueType out[], const otherValueType in[], const IndexType n );

    /** Set out[i] = in[ indexes[i] ],  0 <= i < n */

    template<typename ValueType,typename otherValueType>
    static void setGather( ValueType out[], const otherValueType in[], const IndexType indexes[], const IndexType n );

    /** Set out[ indexes[i] ] = in [i] */

    template<typename ValueType,typename otherValueType>
    static void setScatter( ValueType out[], const IndexType indexes[], const otherValueType in[], const IndexType n );

    /** CUDA implementation for UtilsInterface::Math::invert */

    template<typename ValueType>
    static void invert( ValueType array[], const IndexType n );

    /** Routine that registers all routines of this class at the LAMA interface.
     *
     *  param[inout] UtilsInterface struct to register all routines implemented in CUDA
     */

    static void setInterface( struct UtilsInterface& Utils );

private:

    SCAI_LOG_DECL_STATIC_LOGGER( logger )

    static    bool initialized; //!< static initialization used for registration

    static bool registerInterface();//!< registration
};

} /* end namespace lama */

} /* end namespace scai */
