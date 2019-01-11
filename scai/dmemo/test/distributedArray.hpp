/**
 * @file distributedArray.hpp
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
 * @brief Test Help function to build a distributed array
 * @author Thomas Brandes
 * @date 12.12.2018
 */

namespace scai
{

/** 
 *   Help function to create a distributed array
 *
 *   @param[in] dist is the distribution of the array
 *   @param[in] fill is a function that returns for a global index the value
 *   @returns   the local part owned by this processor of the distributed array
 */
template<typename ValueType>
static hmemo::HArray<ValueType> distributedArray( const dmemo::Distribution& dist, ValueType ( *fill )( IndexType ) )
{
    hmemo::HArray<ValueType> localArray;  // local part of the distributed 'global' array

    // use own scope for write access to make sure that access is closed before return

    {
        IndexType localIndex = 0;   // running local index

        for ( auto& entry : hostWriteOnlyAccess( localArray, dist.getLocalSize() ) )
        {
            entry = fill( dist.local2Global( localIndex++ ) );
        }

    }  // filled the local array with 'global' values

    return localArray;    // each processor gets its local part
}

}
