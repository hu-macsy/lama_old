/**
 * @file ReadAccess.hpp
 *
 * @license
 * Copyright (c) 2011
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
 * @brief Definition of a template class ReadAccess for reading a LAMAArray.
 * @author Thomas Brandes, Jiri Kraus
 * @date 29.04.2011
 * $Id$
 */
#ifndef LAMA_READ_ACCESS_HPP_
#define LAMA_READ_ACCESS_HPP_

// for dll_import
#include <lama/config.hpp>

// base classes
#include <lama/BaseAccess.hpp>

// others
#include <lama/LAMAArray.hpp>
#include <lama/Context.hpp>
#include <lama/LAMAArrayView.hpp>

#include <lama/exception/Exception.hpp>

// logging
#include <logging/logging.hpp>

namespace lama
{

/**
 * @brief The template ReadAccess is used to enforce the consistency of the template LAMAArray.
 *
 * ReadAccess enforces the consistency of the template LAMAArray by following the RAII Idiom. This is
 * done by acquiring a read lock on a LAMAArray in the constructor and releasing this read lock in
 * the destructor. Therefore a ReadAccess should be only used as a stack object.
 */
template<typename T>
class LAMA_DLL_IMPORTEXPORT ReadAccess: public BaseAccess
{
public:
    /**
     * @brief ValueType is the type stored in the wrapped container.
     */
    typedef T ValueType;

    /**
     * @brief acquire a ReadAccess to the passed LAMAArray for the passed context
     *
     * @param[in] array     the LAMAArray to acquire a ReadAccess for
     * @param[in] context   the context that needs a read acess
     * @throws Exception    if the ReadAccess can not be acquired, e.g. because a WriteAccess exists.
     */
    ReadAccess( const LAMAArray<ValueType>& array, ContextPtr context );

    ReadAccess( const LAMAArrayView<ValueType>& view, ContextPtr context );

    ReadAccess( const LAMAArrayConstView<ValueType>& view, ContextPtr context );

    /**
     * @brief Releases the ReadAccess on the associated LAMAArray.
     */
    virtual ~ReadAccess();

    /**
     * @brief returns a valid pointer to the data usable for the context
     *
     * @return a pointer to the wrapped LAMAArray.
     */
    const ValueType* get() const;

    /**
     * @brief release the acquired ReadAccess.
     *
     * Release is mandatory to unlock the array so that it might be
     * used for further write accesses.
     */
    virtual void release();

    virtual void writeAt( std::ostream& stream ) const;

    /**
     * @brief Returns the size of the array
     *
     * @return  the size of the wrapped LAMAArray
     */
    inline IndexType size() const;

private:

    const LAMAArrayConstView<ValueType>* mArrayView;

    size_t mIndex;

    LAMA_LOG_DECL_STATIC_LOGGER(logger);
};

template<typename T>
inline IndexType ReadAccess<T>::size() const
{
    return mArrayView->size();
}

}

#endif // LAMA_READ_ACCESS_HPP_
