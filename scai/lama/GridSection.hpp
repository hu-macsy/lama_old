/**
 * @file GridSector.hpp
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
 * @brief Definition of a multi-dimensional vector.
 * @author Thomas Brandes
 * @date 07.02.2017
 */

#pragma once

#include <scai/logging.hpp>

// for dll_import
#include <scai/common/config.hpp>

#include <scai/common/SCAITypes.hpp>
#include <scai/common/Grid.hpp>
#include <scai/common/BinaryOp.hpp>
#include <scai/common/UnaryOp.hpp>

namespace scai
{

namespace lama
{

// Forward declaration of GridVector required here

template<typename ValueType>
class GridVector;

typedef struct
{
  IndexType globalRange [3];   /* triplet lb:ub:str for global range        */
  IndexType localRange  [3];   /* attention: only valid for not mapped dims */
  bool isRange;                /* false stands for single element           */
  bool isLocalized;            /* true if local range is available          */
} SecDimInfo;

class Range
{
public:

    Range() : mLB( 0 ), mUB( nIndex ), mStride( 1 )
    {
    }

    Range( const IndexType elem ) : mLB( elem ), mUB( elem + 1 ), mStride( 0 )
    {
    }

    Range( const IndexType lb, const IndexType ub ) : mLB( lb ), mUB( ub ), mStride( 1 )
    {
    }

    Range( const IndexType lb, const IndexType ub, const IndexType str ) : mLB( lb ), mUB( ub ), mStride( str )
    {
    }

    IndexType mLB;
    IndexType mUB;
    IndexType mStride;
};

template<typename ValueType>
class COMMON_DLL_IMPORTEXPORT GridSection

{
public:

    GridSection( GridVector<ValueType>& gridVector, const Range& r1 );

    GridSection( GridVector<ValueType>& gridVector, const Range& r1, const Range& r2 );

    GridSection( GridVector<ValueType>& gridVector, const Range& r1, const Range& r2, const Range& r3 );

    GridSection( GridVector<ValueType>& gridVector, const Range& r1, const Range& r2, const Range& r3, const Range& r4 );
     
    template<typename OtherValueType>
    GridSection& operator= ( const GridSection<OtherValueType>& other )
    {
        unaryOp( other, common::unary::COPY );
        return *this;
    }

    GridSection& operator= ( const GridSection& other );

    GridSection& operator*= ( const GridSection& other );

    GridSection& operator+= ( const GridSection& other );

    GridSection& operator-= ( const GridSection& other );

    GridSection& operator/= ( const GridSection& other );

    GridSection& operator= ( const ValueType other );

    GridSection& operator*= ( const ValueType other );

    GridSection& operator+= ( const ValueType other );

    GridSection& operator-= ( const ValueType other );

    GridSection& operator/= ( const ValueType other );

    /** Transpose for two-dimensional sections */

    void assignTranspose( const GridSection<ValueType>& other, bool conjFlag );

    void conj( const GridSection<ValueType>& other );

    void setDim( IndexType k, const Range& r );

public:

    SCAI_LOG_DECL_STATIC_LOGGER( logger )

    GridVector<ValueType>& mGridVector;

    const common::Grid& mGlobalGrid;
    const common::Grid& mLocalGrid;

    const IndexType* localLB;  // lower bounds

    SecDimInfo mDimension[SCAI_GRID_MAX_DIMENSION];
 
    IndexType mNDims;

    bool isConst;  // no write access on section allowed

    /** This operation returns a dope vector of the local part of this section.
     *
     *  @param[out] offset offset for the first element of the section
     *  @param[out] sizes  number of entries for each dimension of the section
     *  @param[out] distances for each dimension
     *  @returns the number of dimensions of this section
     */
    IndexType getDopeVector( IndexType& offset, IndexType sizes[], IndexType distances[] ) const;

    void binOp( const GridSection<ValueType>& other, common::binary::BinaryOp, bool swap );

    void binOpScalar( const ValueType other, common::binary::BinaryOp, bool swap );

    template<typename OtherValueType>
    void unaryOp( const GridSection<OtherValueType>& other, common::unary::UnaryOp op );

    void localize();
};

} /* end namespace lama */

} /* end namespace scai */