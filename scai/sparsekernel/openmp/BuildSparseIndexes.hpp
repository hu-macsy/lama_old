/**
 * @file BuildSparseIndexes.hpp
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
 * @brief Data structure to build sparse indexes
 * @author Thomas Brandes
 * @date 22.12.2016
 */

#pragma once

// for dll_import

#include <scai/common/config.hpp>

#include <scai/common/SCAITypes.hpp>
#include <scai/common/unique_ptr.hpp>

namespace scai
{

namespace sparsekernel
{

/** Data structure that is used to build a sparse vector (or sparse row or col in sparse matrix). 
 *
 *  This data structure is similiar to a stack of indexes.
 *  Indexes already inserted are ignored and this data structure has been optimized for this purpose.
 */

class COMMON_DLL_IMPORTEXPORT BuildSparseIndexes
{
public:

    BuildSparseIndexes( const IndexType n ) : mIndexList( new IndexType[n] )
    {
        NINIT = n + 1; // marks unused colums
        END   = n + 2; // marks end of list

        for ( IndexType j = 0; j < n; j++ )
        {
            mIndexList[j] = NINIT;
        }

        mFirstIndex = END;
        mLength     = 0;
    }

    void pushIndex( const IndexType i )
    {
        if ( mIndexList[i] == NINIT )
        {
            mIndexList[i] = mFirstIndex;
            mFirstIndex = i;
            ++mLength;
        }
    }

    IndexType getLength() const
    {
        return mLength;
    }

    IndexType popIndex() 
    {
        IndexType saveIndex = mFirstIndex;
        
        mFirstIndex = mIndexList[saveIndex];
        mIndexList[saveIndex] = NINIT;

        mLength--;

        return saveIndex;
    }

    bool isEmpty()
    {
        return mFirstIndex == END;
    }

    void reset()
    {
        while ( !isEmpty() )
        {
           popIndex();
        }
    }

private:
 
    common::scoped_array<IndexType> mIndexList;

    IndexType NINIT;
    IndexType END;
    IndexType mFirstIndex;
    IndexType mLength;
};

template<typename ValueType>
class COMMON_DLL_IMPORTEXPORT BuildSparseVector : public BuildSparseIndexes
{
public:

    BuildSparseVector( const IndexType n ) : 

        BuildSparseIndexes( n ),
        mValueList( new ValueType[n] )

    {
        const ValueType zero = 0;

        for ( IndexType j = 0; j < n; j++ )
        {
            mValueList[j] = zero;
        }
    }

    void push( const IndexType i, const ValueType v )
    {
        pushIndex( i );
        mValueList[i] += v;
    }

    void pop( IndexType& i, ValueType& v )
    {
        i = popIndex();
        v = mValueList[i];
        mValueList[i] = static_cast<ValueType>( 0 );
    }

    using BuildSparseIndexes::isEmpty;

private:

    common::scoped_array<ValueType> mValueList;
};


} /* end namespace sparsekernel */

} /* end namespace scai */
