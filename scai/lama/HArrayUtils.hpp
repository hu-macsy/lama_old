/**
 * @file HArrayUtils.hpp
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
 * @brief Definition of class with utility routines.
 * @author Thomas Brandes
 * @date 10.10.2011
 * @since 1.0.0
 */
#pragma once

#include <scai/lama/Scalar.hpp>

// for dll_import
#include <scai/common/config.hpp>

// internal scai libraries
#include <scai/hmemo.hpp>

#include <scai/logging.hpp>

namespace scai
{

namespace lama
{

/** Class that contains some utility routines used at several places. */

class COMMON_DLL_IMPORTEXPORT HArrayUtils
{

public:

    /** Static method for HArray assignment with type conversions.
     *
     *  @param[out] target   contains copy of source values
     *  @param[in]  source   array with source values
     *  @param[in]  context  specifies optionally at which context target will have valid values
     *
     *  \code
     *  HArray<float> fvalues;
     *  HArray<double> dvalues;
     *  ...
     *  dvalues = fvalues;   // not supported
     *  assign( dvalues, fvalues );  // supported
     *  \endcode
     *  Size of target array will be the same as the source array.
     */
    static void assign( hmemo::_HArray& target, const hmemo::_HArray& source, hmemo::ContextPtr context = hmemo::ContextPtr() );

    template<typename ValueType1,typename ValueType2>
    static void assignImpl( hmemo::HArray<ValueType1>& target, const hmemo::HArray<ValueType2>& source, hmemo::ContextPtr context );

    template<typename ValueType1,typename ValueType2>
    static void gather(
        hmemo::HArray<ValueType1>& target,
        const hmemo::HArray<ValueType2>& source,
        const hmemo::HArray<IndexType>& index );

    template<typename ValueType1>
    static void assignScalar( hmemo::HArray<ValueType1>& target, const ValueType1 value, hmemo::ContextPtr context )
                    __attribute__( ( noinline ) );

    static void assignScalar( hmemo::_HArray& target, const Scalar& value, hmemo::ContextPtr context );

    /** This method sets a single value in a heterogeneous array.
     *
     *  @param[in,out] array  Harray where a value to set
     *  @param[in]     index  position to set ( 0 <= index < target.size() )
     *  @param[in]     val    value to set
     *
     *  The value will be set at a valid context.
     */

    template<typename ValueType>
    static void setVal( hmemo::HArray<ValueType>& target, const IndexType index, ValueType val );

    template<typename ValueType>
    static ValueType getVal( const hmemo::HArray<ValueType>& array, const IndexType index );

    /** Scaled assignment on HArray.
     *
     *  @param[out] result  output array
     *  @param[in]  beta    scaling factor
     *  @param[in]  y       source array
     *  @param[in]  context location where operation is done
     */

    template<typename ValueType>
    static void assignScaled(
        hmemo::HArray<ValueType>& result,
        const ValueType beta,
        const hmemo::HArray<ValueType>& y,
        hmemo::ContextPtr context );

    /** scale array in place 
     *
     *  Note: scale will be done where array has currently valid values. The preferred
     *        location is not taken if the array is not valid there.
     */

    template<typename ValueType>
    static void scale( hmemo::HArray<ValueType>& array, const ValueType beta, hmemo::ContextPtr prefLoc );

    /** Replace in a complex array its values with the conjugate values */

    template<typename ValueType>
    static void conj( hmemo::HArray<ValueType>& array, hmemo::ContextPtr prefLoc );

private:

    template<typename ValueType>
    static void assignImpl1( hmemo::HArray<ValueType>& target, const hmemo::_HArray& source, hmemo::ContextPtr context );

    SCAI_LOG_DECL_STATIC_LOGGER( logger )
};

} /* end namespace lama */

} /* end namespace scai */
