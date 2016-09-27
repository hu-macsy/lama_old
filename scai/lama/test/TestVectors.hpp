/**
 * @file TestVectors.hpp
 *
 * @license
 * Copyright (c) 2009-2016
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
 * @brief Array with all vectors, one for each supported vector storage kind/type
 * @author Thomas Brandes
 * @date 27.07.2016
 */

#include <scai/lama/Vector.hpp>

#include <scai/common/shared_ptr.hpp>
#include <scai/common/TypeTraits.hpp>

#include <vector>

namespace scai
{

namespace lama
{

/** Class for a list of vectors, one for each supported
 *  vector storage format and each supported arithmetic type.
 *  
 *  Note: Currently only DENSE vectors are supported, but SPARSE might be 
 *        supported in future releases.
 *  
 */

class TestVectors : public std::vector<scai::lama::VectorPtr>
{

public:

    /** Constructor creates already the list with shared vector pointers, one for each
     *  registered key type in the Vector factory.
     *
     *  @param[in] ctx optional argument for the context where operations on vector should be executed
     */

    TestVectors( scai::hmemo::ContextPtr ctx = scai::hmemo::ContextPtr() )
    {
        using namespace scai::lama;
        std::vector<VectorCreateKeyType> values;  //  all create values
        Vector::getCreateValues( values );

        for ( size_t i = 0; i < values.size(); ++i )
        {
            scai::lama::VectorPtr vectorPtr( scai::lama::Vector::create( values[i] ) );

            if ( ctx )
            {
                vectorPtr->setContextPtr( ctx );
            }

            push_back( vectorPtr );
        }
    }
};

}

}

