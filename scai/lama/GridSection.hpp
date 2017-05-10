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

// for dll_import
#include <scai/common/config.hpp>


#include <scai/lama/GridVector.hpp>

namespace scai
{

namespace lama
{

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

    GridSection( GridVector<ValueType>& gridVector, const Range& r1, const Range& r2, const Range& r3 ) :
     
        mGridVector( gridVector ), 
        mGlobalGrid( gridVector.globalGrid() ),
        mLocalGrid( gridVector.localGrid() )

    {
        mNDims = 3;

        setDim( 0, r1 );
        setDim( 1, r2 );
        setDim( 2, r3 );
    }

    GridSection& operator= ( const GridSection& other );

    GridSection& operator= ( const ValueType other );

    GridSection& operator*= ( const ValueType other );

    void setDim( IndexType k, const Range& r )
    {
       SecDimInfo& sec = mDimension[k];

       sec.globalRange[0] = r.mLB;
       sec.globalRange[1] = r.mUB;
       sec.globalRange[2] = r.mStride;

       sec.isLocalized = false;
       sec.isRange     = r.mStride != 0;

       if ( r.mUB == nIndex )
       {
           sec.globalRange[1] = mGlobalGrid.size( k );
       }

       if ( r.mStride != 0 )
       {
           sec.isRange = true;
           sec.globalRange[2] = r.mStride;
       }
       else
       {
           sec.isRange = false;
           sec.globalRange[2] = 1;
       }
    }

private:

    GridVector<ValueType>& mGridVector;

    const common::Grid& mGlobalGrid;
    const common::Grid& mLocalGrid;

    SecDimInfo mDimension[SCAI_GRID_MAX_DIMENSION];
 
    IndexType mNDims;

    bool isConst;  // no write access on section allowed
};

} /* end namespace lama */

} /* end namespace scai */
