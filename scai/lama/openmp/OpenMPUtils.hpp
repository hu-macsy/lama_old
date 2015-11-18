/**
 * @file OpenMPUtils.hpp
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
 * @brief Implementation of general utilities with OpenMP
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
#include <scai/common/macros/assert.hpp>

namespace scai
{

namespace lama
{

/** General utilities of the LAMA Interface implemented in OpenMP  */

class COMMON_DLL_IMPORTEXPORT OpenMPUtils
{
public:

    /** OpenMP implementation for UtilKernelTrait::Transform::scale */

    template<typename ValueType>
    static void scale( ValueType mValues[], const ValueType value, const IndexType n );

    /** OpenMP implementation for UtilKernelTrait::Copy::setScale */

    template<typename ValueType,typename OtherValueType>
    static void setScale(
        ValueType outValues[],
        const ValueType value,
        const OtherValueType inValues[],
        const IndexType n );

    /*  This method is an implementation of UtilKernelTrait::validIndexes */

    static bool validIndexes( const IndexType array[], const IndexType n, const IndexType size );

    /** OpenMP implementation for UtilKernelTrait::Reductions::sum */

    template<typename ValueType>
    static ValueType sum( const ValueType array[], const IndexType n );

    /** OpenMP implementation for UtilKernelTrait::Setter::setVal */

    template<typename ValueType>
    static void setVal( ValueType array[], const IndexType n, const ValueType val );

    /** OpenMP implementation for UtilKernelTrait::Setter::setOrder */

    template<typename ValueType>
    static void setOrder( ValueType array[], const IndexType n );

    template<typename ValueType>
    static ValueType getValue( const ValueType* array, const IndexType i );

    template<typename ValueType>
    static ValueType maxval( const ValueType array[], const IndexType n );

    /** OpenMP implementation for UtilKernelTrait::Reductions::absMaxVal */

    template<typename ValueType>
    static ValueType absMaxVal( const ValueType array[], const IndexType n );

    /** OpenMP implementation for UtilKernelTrait::Reductions::absMaxDiffVal */

    template<typename ValueType>
    static ValueType absMaxDiffVal( const ValueType array1[], const ValueType array2[], const IndexType n );

    /** OpenMP implementation for UtilKernelTrait::Reductions::isSorted */

    template<typename ValueType>
    static bool isSorted( const ValueType array[], const IndexType n, bool acending );

    template<typename ValueType1,typename ValueType2>
    static void set( ValueType1 out[], const ValueType2 in[], const IndexType n );

    /** Set out[i] = in[ indexes[i] ],  0 <= i < n */

    template<typename ValueType1,typename ValueType2>
    static void setGather( ValueType1 out[], const ValueType2 in[], const IndexType indexes[], const IndexType n );

    /** Set out[ indexes[i] ] = in [i] */

    template<typename ValueType1,typename ValueType2>
    static void setScatter( ValueType1 out[], const IndexType indexes[], const ValueType2 in[], const IndexType n );

    /** OpenMP implementation for UtilKernelTrait::Math::invert */

    template<typename ValueType>
    static void invert( ValueType array[], const IndexType n );

private:

    /** Routine that registers all methods at the kernel registry. */

    static void registerKernels( bool deleteFlag );

    /** Constructor for registration. */

    OpenMPUtils();

    /** Destructor for unregistration. */

    ~OpenMPUtils();

    /** Static variable for registration at static initialization. */

    static OpenMPUtils guard;

    SCAI_LOG_DECL_STATIC_LOGGER( logger )
};

} /* end namespace lama */

} /* end namespace scai */
