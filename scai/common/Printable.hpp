/**
 * @file Printable.hpp
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
 * @brief Base class to be used for base classes that will output on a stream.
 * @author Jiri Kraus
 * @date 05.08.2011
 */

#pragma once

// for dll_import
#include <scai/common/config.hpp>

#include <ostream>

namespace scai
{

namespace common
{

/**
 * @brief The class Printable is used as a base class for all classes
 *        that should provide the output operator << to print info about it.
 *
 * This base class avoids the definition of the operator<< for individual
 * classes.
 *
 * The operator<< is especially helpful for logging messages.
 */
class COMMON_DLL_IMPORTEXPORT Printable
{
public:
    /**
     * @brief Creates a Printable.
     */
    Printable();

    /**
     * @brief Destroys a Printable.
     */
    virtual ~Printable();

    /**
     * @brief Writes some information about this to the passed stream.
     *
     * The method should be overwritten by base classes to give more
     * specific information about the object.
     * If a deriving class does not override it, typeid(this).name() is
     * written to stream.
     *
     * @param[out]  stream  the stream to write to.
     */
    virtual void writeAt( std::ostream& stream ) const;

    /**
     * @brief Sets or resets extended output mode for printing structures.
     */
    static void enableExtended( const bool flag );

protected:

    static bool extended; //!< flag for extended output of objects.
};

/**
 * @brief Calls Printable::writeAt() on the passed stream.
 *
 * @param[out] stream   the stream to write to.
 * @param[in]  object   the Printable which should write to stream.
 *
 * @return a reference to the passed stream.
 */

COMMON_DLL_IMPORTEXPORT std::ostream& operator<<( std::ostream& stream, const Printable& object );

} /* end namespace common */

} /* end namespace scai */
