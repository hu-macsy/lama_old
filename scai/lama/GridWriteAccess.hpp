/**
 * @file GridWriteAccess.hpp
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
 * @brief Derived class of WriteAccess to write local data of a GridVector
 * @author Thomas Brandes
 * @date 10.05.2017
 */

#pragma once

// for dll_import
#include <scai/common/config.hpp>

// base classes
#include <scai/hmemo/WriteAccess.hpp>
#include <scai/lama/GridVector.hpp>

namespace scai
{

namespace lama
{

/**
 * @brief Write access for the local data of a grid vector.
 *
 *
 * @tparam ValueType is the value type of the data written.
 */
template<typename ValueType>
class COMMON_DLL_IMPORTEXPORT GridWriteAccess: public hmemo::WriteAccess<ValueType>
{
public:

    GridWriteAccess( GridVector<ValueType>& gridVector, hmemo::ContextPtr contextPtr );

    GridWriteAccess( GridVector<ValueType>& gridVector );

    /**
     * @brief Releases the WriteAccess on the associated GridVector
     */
    virtual ~GridWriteAccess();

    /**
     *  Indexing a one-dimensional grid is like a usual write access but added for convenience.
     */
    ValueType& operator() ( const IndexType i1 );

    /**
     *  Indexing a two-dimensional grid encapsulates the linear addressing that is required.
     */
    ValueType& operator() ( const IndexType i1, const IndexType i2 );

    ValueType& operator() ( const IndexType i1, const IndexType i2, const IndexType i3 );

    ValueType& operator() ( const IndexType i1, const IndexType i2, const IndexType i3, const IndexType i4 );

    IndexType size( const IndexType dim );

    using hmemo::WriteAccess<ValueType>::size;    // make size() visible

private:
 
    const common::Grid& mGrid;
};

/* --------------------------------------------------------------------------- */

template<typename ValueType>
GridWriteAccess<ValueType>::GridWriteAccess( GridVector<ValueType>& gv, hmemo::ContextPtr contextPtr ) :

    hmemo::WriteAccess<ValueType>( gv.getLocalValues(), contextPtr ),
    mGrid( gv.localGrid() )
{
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
GridWriteAccess<ValueType>::GridWriteAccess( GridVector<ValueType>& gv ) :

    hmemo::WriteAccess<ValueType>( gv.getLocalValues() ),
    mGrid( gv.localGrid() )
{
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
GridWriteAccess<ValueType>::~GridWriteAccess()
{
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
IndexType GridWriteAccess<ValueType>::size( const IndexType dim )
{
    return mGrid.size( dim );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
ValueType& GridWriteAccess<ValueType>::operator() ( const IndexType i1 )
{
    SCAI_ASSERT_EQ_DEBUG( 1, mGrid.nDims(), "illegal indexing" )
    SCAI_ASSERT_VALID_INDEX_DEBUG( i1, this->size(), "out of range pos" )
    return (*this)[i1];
}

template<typename ValueType>
ValueType& GridWriteAccess<ValueType>::operator() ( const IndexType i1, const IndexType i2 )
{
    SCAI_ASSERT_EQ_DEBUG( 2, mGrid.nDims(), "illegal indexing" )
    const common::Grid2D& grid = reinterpret_cast<const common::Grid2D&>( mGrid );
    IndexType pos = grid.linearPos( i1, i2 );
    SCAI_ASSERT_VALID_INDEX_DEBUG( pos, this->size(), "out of range pos" )
    return (*this)[pos];
}

template<typename ValueType>
ValueType& GridWriteAccess<ValueType>::operator() ( const IndexType i1, const IndexType i2, const IndexType i3 )
{
    SCAI_ASSERT_EQ_DEBUG( 3, mGrid.nDims(), "illegal indexing" )
    const common::Grid3D& grid = reinterpret_cast<const common::Grid3D&>( mGrid );
    IndexType pos = grid.linearPos( i1, i2, i3 );
    SCAI_ASSERT_VALID_INDEX_DEBUG( pos, this->size(), "out of range pos" )
    return (*this)[pos];
}

template<typename ValueType>
ValueType& GridWriteAccess<ValueType>::operator() ( const IndexType i1, const IndexType i2, const IndexType i3, IndexType i4 )
{
    SCAI_ASSERT_EQ_DEBUG( 4, mGrid.nDims(), "illegal indexing" )
    const common::Grid4D& grid = reinterpret_cast<const common::Grid4D&>( mGrid );
    IndexType pos = grid.linearPos( i1, i2, i3, i4 );
    SCAI_ASSERT_VALID_INDEX_DEBUG( pos, this->size(), "out of range pos (" << i1 << ", " << i2 << ", " << i3 << ", " << i4 << " of " << grid )
    return (*this)[pos];
}

} /* end namespace lama */

} /* end namespace scai */
