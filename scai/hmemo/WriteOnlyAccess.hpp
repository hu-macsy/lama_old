/**
 * @file WriteOnlyAccess.hpp
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
 * @brief Definition of a template class for rewriting of a LAMA array.
 * @author Thomas Brandes
 * @date 16.07.2015
 */

#pragma once

// for dll_import
#include <scai/common/config.hpp>

// base classes
#include <scai/hmemo/WriteAccess.hpp>

// logging
#include <scai/logging.hpp>
#include <scai/common/TypeTraits.hpp>

namespace scai
{

namespace hmemo
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

    WriteOnlyAccess( HArray<ValueType>& array, ContextPtr context )
        : WriteAccess<ValueType>( array, context, false )
    {
        SCAI_LOG_DEBUG( logger, "WriteOnlyAccess<" << common::TypeTraits<ValueType>::id() << ">" )
    }

    /** Create a write access with keep flag = false for the host context
     *
     *  Attention: this kind of write access assumes that the array is completely new written.
     */

    WriteOnlyAccess( HArray<ValueType>& array )
        : WriteAccess<ValueType>( array, false )
    {
        SCAI_LOG_DEBUG( logger, "WriteOnlyAccess<" << common::TypeTraits<ValueType>::id() << ">" )
    }

    /**
     *  @brief Acquire a write only access with a resize.
     *
     * @param[in] array     the HArray for which access is required
     * @param[in] context   the context where data will be written later
     * @param[in] size      the new size of the LAMA array.
     *
     * Attention: this kind of write access assumes that the array is completely new written.
     */
    WriteOnlyAccess( HArray<ValueType>& array, ContextPtr context, const IndexType size )
        : WriteAccess<ValueType>( array, context, false )
    {
        this->resize( 0 );      // invalidates all data before resize
        this->resize( size );   // now resize
        SCAI_LOG_DEBUG( logger, "WriteOnlyAccess<" << common::TypeTraits<ValueType>::id() << ">: " << *mArray )
    }

    WriteOnlyAccess( HArray<ValueType>& array, const IndexType size )
        : WriteAccess<ValueType>( array, false )
    {
        this->resize( 0 );      // invalidates all data before resize
        this->resize( size );   // now resize
        SCAI_LOG_DEBUG( logger, "WriteOnlyAccess<" << common::TypeTraits<ValueType>::id() << ">: " << *mArray )
    }

    WriteOnlyAccess ( WriteOnlyAccess<ValueType>&& other )
        : WriteAccess<ValueType>( std::move( other ) )
    { }

    ~WriteOnlyAccess()
    {
        SCAI_LOG_DEBUG( WriteAccess<ValueType>::logger, "~WriteOnlyAccess<" << common::TypeTraits<ValueType>::id() << ">" )
    }

protected:

#ifndef SCAI_LOG_LEVEL_OFF
    using WriteAccess<ValueType>::logger;   // no additinal logger for this derived class
#endif

private:

    using WriteAccess<ValueType>::mArray;
    using WriteAccess<ValueType>::mContextDataIndex;
};

/**
 * @brief Obtain a WriteOnlyAccess to the given array.
 *
 * Analogous to readAccess(const HArray & array). See its documentation for
 * motivation and intended usage of this function.
 */
template <typename ValueType>
WriteOnlyAccess<ValueType> writeOnlyAccess( HArray<ValueType>& array, IndexType newSize )
{
    return WriteOnlyAccess<ValueType>( array, newSize );
}

/**
 * @brief Obtain a WriteOnlyAccess for the supplied context to the given array.
 *
 * Analogous to readAccess(const HArray & array, ContextPtr). See its documentation for
 * motivation and intended usage of this function.
 */
template <typename ValueType>
WriteOnlyAccess<ValueType> writeOnlyAccess( HArray<ValueType>& array, ContextPtr context, IndexType newSize )
{
    return WriteOnlyAccess<ValueType>( array, context, newSize );
}

} /* end namespace hmemo */

} /* end namespace scai */
