/**
 * @file LAMAArrayUtils.hpp
 *
 * @license
 * Copyright (c) 2009-2013
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
#ifndef LAMA_LAMA_ARRAY_UTILS_HPP_
#define LAMA_LAMA_ARRAY_UTILS_HPP_
// for dll_import
#include <lama/config.hpp>

// others
#include <lama/LAMAArray.hpp>

// logging
#include <logging/logging.hpp>

namespace lama
{

/** Class that contains some utility routines used at several places. */

class LAMA_DLL_IMPORTEXPORT LAMAArrayUtils
{

public:

    /** Static method for LAMAArray assignment with type conversions.
     *
     *  @param[out] target   contains copy of source values
     *  @param[in]  source   array with source values
     *  @param[in]  context  specifies optionally at which context target will have valid values
     *
     *  \code
     *  LAMAArray<float> fvalues;
     *  LAMAArray<double> dvalues;
     *  ...
     *  dvalues = fvalues;   // not supported
     *  assign( dvalues, fvalues );  // supported
     *  \endcode
     *  Size of target array will be the same as the source array.
     */
    static void assign( _LAMAArray& target, const _LAMAArray& source, ContextPtr context = ContextPtr() );

    template<typename ValueType1,typename ValueType2>
    static void assignImpl( LAMAArray<ValueType1>& target, const LAMAArray<ValueType2>& source, ContextPtr context );

    template<typename ValueType1,typename ValueType2>
    static void gather(
        LAMAArray<ValueType1>& target,
        const LAMAArray<ValueType2>& source,
        const LAMAArray<IndexType>& index );

    template<typename ValueType1>
    static void assignScalar( LAMAArray<ValueType1>& target, const Scalar& value, ContextPtr context );

    static void assignScalar( _LAMAArray& target, const Scalar& value, ContextPtr context );

    /** This method sets a single value in a LAMA array.
     *
     *  @param[in,out] target LAMA array where a value to set
     *  @param[in]     index  position to set ( 0 <= index < target.size() )
     *  @param[in]     val    value to set
     */

    template<typename ValueType>
    static void setVal( LAMAArray<ValueType>& target, const IndexType index, ValueType val );

private:

    template<typename ValueType>
    static void assignImpl1( LAMAArray<ValueType>& target, const _LAMAArray& source, ContextPtr context );

    LAMA_LOG_DECL_STATIC_LOGGER( logger )
};

} // namespace

#endif // LAMA_LAMA_ARRAY_UTILS_HPP_
