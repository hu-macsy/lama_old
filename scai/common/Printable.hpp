/**
 * @file Printable.hpp
 *
 * @license
 * Copyright (c) 2009-2017
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the SCAI framework LAMA.
 *
 * LAMA is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Affero General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * LAMA is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for
 * more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with LAMA. If not, see <http://www.gnu.org/licenses/>.
 *
 * Other Usage
 * Alternatively, this file may be used in accordance with the terms and
 * conditions contained in a signed written agreement between you and
 * Fraunhofer SCAI. Please contact our distributor via info[at]scapos.com.
 * @endlicense
 *
 * @brief Base class to be used for base classes that will output on a stream.
 * @author Jiri Kraus
 * @date 05.08.2011
 */

#pragma once

// for dll_import
#include <scai/common/config.hpp>

// std
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

/**
 *  @brief Wrapper class to print contents of container
 *
 *  \code
 *     std::vector<int> list( { 0, 1, 2, 3, 5, 7 } );
 *     std::cout << "Content of list = " << common::Wrapper( list ) << std::endl;
 *  \endcode
 */
template<class Container>
class Wrapper : public Printable
{

public:

    Wrapper( const Container& container ) :  mContainer( container ) 
    {
    }

    ~Wrapper()
    {
    }

    virtual void writeAt( std::ostream& stream ) const
    {
        typename Container::const_iterator beg = mContainer.begin();

        stream << "["; 

        while( beg != mContainer.end())
        {
            stream << " " << *beg++; 
        }

        stream  << " ]"; 
    }

private:

    const Container& mContainer;  // keep a reference
};

} /* end namespace common */

} /* end namespace scai */
