/**
 * @file GridReadAccess.hpp
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
 * @brief Derived class of ReadAccess to write local data of a GridVector
 * @author Thomas Brandes
 * @date 10.05.2017
 */

#pragma once

// for dll_import
#include <scai/common/config.hpp>

// base classes
#include <scai/hmemo/ReadAccess.hpp>
#include <scai/lama/GridVector.hpp>

namespace scai
{

namespace lama
{

/**
 * @brief Read access for the local data of a grid vector.
 *
 * This class derives from ReadAccess but keeps also the information about
 * the local grid sizes. By this way it allows multi-dimensinal indexing
 * that is less error prone than linearized indexes. 
 *
 * @tparam ValueType is the value type of the data written.
 */
template<typename ValueType>
class COMMON_DLL_IMPORTEXPORT GridReadAccess: public hmemo::ReadAccess<ValueType>
{
public:

    /** Constructor by a grid vector and the desired context where the local grid data is read. */

    GridReadAccess( const GridVector<ValueType>& gridVector, hmemo::ContextPtr contextPtr );

    /** Constructor by a grid vector, the access defaults to the Host. */

    GridReadAccess( const GridVector<ValueType>& gridVector );

    /**
     * @brief Releases the ReadAccess on the associated GridVector
     */
    virtual ~GridReadAccess();

    const ValueType& operator() ( const IndexType i1 );

    const ValueType& operator() ( const IndexType i1, const IndexType i2 );

    const ValueType& operator() ( const IndexType i1, const IndexType i2, const IndexType i3 );

    const ValueType& operator() ( const IndexType i1, const IndexType i2, const IndexType i3, const IndexType i4 );

private:
 
    const common::Grid& mGrid;  // keep a reference to grid
};

/* --------------------------------------------------------------------------- */

template<typename ValueType>
GridReadAccess<ValueType>::GridReadAccess( const GridVector<ValueType>& gv, hmemo::ContextPtr contextPtr ) :

    hmemo::ReadAccess<ValueType>( gv.getLocalValues(), contextPtr ),
    mGrid( gv.localGrid() )
{
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
GridReadAccess<ValueType>::GridReadAccess( const GridVector<ValueType>& gv ) :

    hmemo::ReadAccess<ValueType>( gv.getLocalValues() ),
    mGrid( gv.localGrid() )
{
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
GridReadAccess<ValueType>::~GridReadAccess()
{
}

template<typename ValueType>
const ValueType& GridReadAccess<ValueType>::operator() ( const IndexType i1 )
{
    SCAI_ASSERT_EQ_DEBUG( 1, mGrid.nDims(), "illegal indexing" )
    // call of linear pos redundant, is just identity
    return (*this)[i1];
}

template<typename ValueType>
const ValueType& GridReadAccess<ValueType>::operator() ( const IndexType i1, const IndexType i2 )
{
    SCAI_ASSERT_EQ_DEBUG( 2, mGrid.nDims(), "illegal indexing" )
    const common::Grid2D& grid = reinterpret_cast<const common::Grid2D&>( mGrid );
    IndexType pos = grid.linearPos( i1, i2 );
    return (*this)[pos];
}

template<typename ValueType>
const ValueType& GridReadAccess<ValueType>::operator() ( const IndexType i1, const IndexType i2, const IndexType i3 )
{
    SCAI_ASSERT_EQ_DEBUG( 3, mGrid.nDims(), "illegal indexing" )
    const common::Grid3D& grid = reinterpret_cast<const common::Grid3D&>( mGrid );
    IndexType pos = grid.linearPos( i1, i2, i3 );
    return (*this)[pos];
}

template<typename ValueType>
const ValueType& GridReadAccess<ValueType>::operator() ( const IndexType i1, const IndexType i2, const IndexType i3, const IndexType i4 )
{
    SCAI_ASSERT_EQ_DEBUG( 4, mGrid.nDims(), "illegal indexing" )
    const common::Grid4D& grid = reinterpret_cast<const common::Grid4D&>( mGrid );
    IndexType pos = grid.linearPos( i1, i2, i3, i4 );
    return (*this)[pos];
}

} /* end namespace lama */

} /* end namespace scai */
