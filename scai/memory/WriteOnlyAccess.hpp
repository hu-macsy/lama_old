/**
 * @file WriteOnlyAccess.hpp
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
 * @brief Definition of a template class for rewriting of a LAMA array.
 * @author Thomas Brandes
 * @date 16.07.2015
 */

#pragma once

// for dll_import
#include <scai/common/config.hpp>

// base classes
#include <scai/memory/WriteAccess.hpp>

// logging
#include <scai/logging.hpp>

namespace scai
{

namespace memory
{

/**
 * @brief WriteOnlyAccess is a write access where no existing values of the array are needed (keepFlag).
 *
 * This derived class has been added for more convenience as it avoids the use of the keepFlag param.
 *
 * A WriteOnlyAccess should be used whenever possible. It avoids any memory transfer of no more
 * needed values between devices and in case of a reallocation it avoids copying of old values.
 *
 * @tparam ValueType is the value type stored in the wrapped container.
 */
template<typename ValueType>
class COMMON_DLL_IMPORTEXPORT WriteOnlyAccess: public WriteAccess<ValueType>
{
public:

    /** Create a write access with keep flag = false.
     *
     *  Attention: this kind of write access assumes that the array is completely new written.
     */

    WriteOnlyAccess( LAMAArray<ValueType>& array, ContextPtr context )
        : WriteAccess<ValueType>( array, context, false )
    {
        SCAI_LOG_DEBUG( logger, "WriteOnlyAccess<" << common::getScalarType<ValueType>() << ">" )
    }

    /** Create a write access with keep flag = false for the host context
     *
     *  Attention: this kind of write access assumes that the array is completely new written.
     */

    WriteOnlyAccess( LAMAArray<ValueType>& array )
        : WriteAccess<ValueType>( array, false )
    {
        SCAI_LOG_DEBUG( logger, "WriteOnlyAccess<" << common::getScalarType<ValueType>() << ">" )
    }

    /**
     *  @brief Acquire a write only access with a resize.
     *
     * @param[in] array     the LAMAArray for which access is required
     * @param[in] context   the context where data will be written later
     * @param[in] size      the new size of the LAMA array.
     *
     * Attention: this kind of write access assumes that the array is completely new written.
     */
    WriteOnlyAccess( LAMAArray<ValueType>& array, ContextPtr context, const common::IndexType size )
        : WriteAccess<ValueType>( array, context, false )
    {
        this->resize( 0 );      // invalidates all data before resize
        this->resize( size );   // now resize

        SCAI_LOG_DEBUG( logger, "WriteOnlyAccess<" << common::getScalarType<ValueType>() << ">: " << *mArray )
    }

    WriteOnlyAccess( LAMAArray<ValueType>& array, const common::IndexType size )
        : WriteAccess<ValueType>( array, false )
    {
        this->resize( 0 );      // invalidates all data before resize
        this->resize( size );   // now resize

        SCAI_LOG_DEBUG( logger, "WriteOnlyAccess<" << common::getScalarType<ValueType>() << ">: " << *mArray )
    }

    ~WriteOnlyAccess()
    {
        SCAI_LOG_DEBUG( WriteAccess<ValueType>::logger, "~WriteOnlyAccess<" << common::getScalarType<ValueType>() << ">" )
    }

protected:

    using WriteAccess<ValueType>::logger;   // no additinal logger for this derived class

private:

    using WriteAccess<ValueType>::mArray;
    using WriteAccess<ValueType>::mContextDataIndex;
};

} /* end namespace memory */

} /* end namespace scai */
