/**
 * @file HArrays.hpp
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
 * @brief Vector with all HArray one for each supported type
 * @author Thomas Brandes
 * @date 31.08.2012
 */

#include <scai/hmemo/_HArray.hpp>

#include <vector>

typedef std::shared_ptr<scai::hmemo::_HArray> ArrayPtr;

/** Class for a list of arrays, one for each supported array type.
 */

class HArrays : public std::vector<ArrayPtr>
{

public:

    /** Constructor creates already the list with all storage pointers. */

    HArrays( scai::hmemo::ContextPtr ctx = scai::hmemo::ContextPtr() )
    {
        using namespace scai::common;
        using namespace scai::hmemo;
        std::vector<ScalarType> values;  //  all create values
        _HArray::getCreateValues( values );

        for ( size_t i = 0; i < values.size(); ++i )
        {
            ArrayPtr arrayPtr( _HArray::create( values[i] ) );

            if ( ctx )
            {
                arrayPtr->prefetch( ctx );
            }

            push_back( arrayPtr );
        }
    }

    // Destructor will free all arrays due to use of shared pointers
};
