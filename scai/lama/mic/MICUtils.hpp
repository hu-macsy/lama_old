/**
 * @file MICUtils.hpp
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
 * @brief Implementation of general utilities with MIC
 * @author Thomas Brandes
 * @date 02.07.2012
 * @since 1.1.0
 */

#pragma once

// for dll_import
#include <scai/common/config.hpp>

// others
#include <scai/common/SCAITypes.hpp>
#include <scai/common/Assert.hpp>
#include <scai/hmemo/mic/MICContext.hpp>

// logging
#include <scai/logging.hpp>

namespace scai
{

using common::IndexType;

namespace lama
{

/** General utilities of the LAMA Interface implemented in MIC  */

class COMMON_DLL_IMPORTEXPORT MICUtils
{
public:

    /** MIC implementation for UtilsInterface::Transform::scale */

    template<typename ValueType>
    static void scale( ValueType mValues[], const ValueType value, const IndexType n );

    /** MIC implementation for UtilsInterface::Copy::setScale */

    template<typename ValueType,typename OtherValueType>
    static void setScale(
        ValueType outValues[],
        const ValueType value,
        const OtherValueType inValues[],
        const IndexType n );

    /*  This method is an implementation of UtilsInterface::validIndexes */

    static bool validIndexes( const IndexType array[], const IndexType n, const IndexType size );

    /** MIC implementation for UtilsInterface::Reductions::sum */

    template<typename ValueType>
    static ValueType sum( const ValueType array[], const IndexType n );

    /** MIC implementation for UtilsInterface::Setter::setVal */

    template<typename ValueType>
    static void setVal( ValueType array[], const IndexType n, const ValueType val );

    /** MIC implementation for UtilsInterface::Setter::setOrder */

    template<typename ValueType>
    static void setOrder( ValueType array[], const IndexType n );

    template<typename ValueType>
    static ValueType getValue( const ValueType* array, const IndexType i );

    template<typename ValueType>
    static ValueType maxval( const ValueType array[], const IndexType n );

    /** MIC implementation for UtilsInterface::Reductions::absMaxVal */

    template<typename ValueType>
    static ValueType absMaxVal( const ValueType array[], const IndexType n );

    /** MIC implementation for UtilsInterface::Reductions::absMaxDiffVal */

    template<typename ValueType>
    static ValueType absMaxDiffVal( const ValueType array1[], const ValueType array2[], const IndexType n );

    /** MIC implementation for UtilsInterface::Reductions::isSorted */

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

    /** MIC implementation for UtilsInterface::Math::invert */

    template<typename ValueType>
    static void invert( ValueType array[], const IndexType n );

private:

    /** Routine that registers all routines of this class at the LAMA interface.
     *
     *  param[inout] UtilsInterface struct to register all routines implemented in MIC
     */

    static void setInterface( struct UtilsInterface& Utils );

    static bool initialized;

    static bool registerInterface();

    SCAI_LOG_DECL_STATIC_LOGGER( logger )
};

} /* end namespace lama */

} /* end namespace scai */
