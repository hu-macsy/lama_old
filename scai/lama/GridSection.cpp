/**
 * @file GridSection.cpp
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
 * @brief Implementation of operations on distributed sections
 * @author Thomas Brandes
 * @date 10.05.2017
 */

#include <scai/lama/GridSection.hpp>
#include <scai/lama/GridVector.hpp>

#include <scai/lama/GridReadAccess.hpp>
#include <scai/lama/GridWriteAccess.hpp>

#include <scai/common/macros/instantiate.hpp>
#include <scai/common/SCAITypes.hpp>

namespace scai
{

namespace lama
{

/* ---------------------------------------------------------------------------------------*/

SCAI_LOG_DEF_TEMPLATE_LOGGER( template<typename ValueType>, GridSection<ValueType>::logger, "GridSection" )

/* ---------------------------------------------------------------------------------------*/

template<typename ValueType>
GridSection<ValueType>::GridSection( GridVector<ValueType>& gridVector, const Range& r0 ) :

    mGridVector( gridVector ),
    mGlobalGrid( gridVector.globalGrid() ),
    mLocalGrid( gridVector.localGrid() )

{
    mNDims = 1;

    setDim( 0, r0 );
}

template<typename ValueType>
GridSection<ValueType>::GridSection( GridVector<ValueType>& gridVector, const Range& r0, const Range& r1 ) :

    mGridVector( gridVector ),
    mGlobalGrid( gridVector.globalGrid() ),
    mLocalGrid( gridVector.localGrid() )
{
    mNDims = 2;

    setDim( 0, r0 );
    setDim( 1, r1 );
}

template<typename ValueType>
GridSection<ValueType>::GridSection( GridVector<ValueType>& gridVector, const Range& r0, const Range& r1, const Range& r2 ) :

    mGridVector( gridVector ),
    mGlobalGrid( gridVector.globalGrid() ),
    mLocalGrid( gridVector.localGrid() )
{
    mNDims = 3;

    setDim( 0, r0 );
    setDim( 1, r1 );
    setDim( 2, r2 );
}

template<typename ValueType>
GridSection<ValueType>::GridSection( GridVector<ValueType>& gridVector, const Range& r0, const Range& r1, const Range& r2, const Range& r3 ) :

    mGridVector( gridVector ),
    mGlobalGrid( gridVector.globalGrid() ),
    mLocalGrid( gridVector.localGrid() )
{
    mNDims = 4;

    setDim( 0, r0 );
    setDim( 1, r1 );
    setDim( 2, r2 );
    setDim( 3, r3 );
}

template<typename ValueType>
IndexType GridSection<ValueType>::getDopeVector( IndexType& offset, IndexType sizes[], IndexType distances[] ) const
{
    SCAI_LOG_ERROR( logger, "get dope vector grid section, #dims = " << mNDims )

    mLocalGrid.getDistances( distances );

    offset = 0;
    
    IndexType nSectionDims = 0;  // counts only ranges, not fixed elements

    for ( IndexType i = 0; i < mNDims; ++i )
    {
        const SecDimInfo& sec = mDimension[i];

        SCAI_LOG_ERROR( logger, "secDim[" << i << "] : " << sec.globalRange[0] << ":" << sec.globalRange[1] << ":" << sec.globalRange[2] 
                                 << ", is range = " << sec.isRange )

        offset += sec.globalRange[0] * distances[i];

        SCAI_LOG_ERROR( logger, "offset = " << offset )

        if ( sec.isRange )
        {
            SCAI_ASSERT_EQ_ERROR( 1, sec.globalRange[2], "stride != 1 not supported yet" )
            sizes[ nSectionDims ] = sec.globalRange[1] - sec.globalRange[0];
            distances[ nSectionDims] = distances[ i ];
            nSectionDims++;
        }
    }

    return nSectionDims;
}

template<typename ValueType>
GridSection<ValueType>& GridSection<ValueType>::operator= ( const GridSection<ValueType>& other )
{
    IndexType offsetSource = 0;
    IndexType offsetTarget = 0;

    IndexType sizesSource[SCAI_GRID_MAX_DIMENSION];
    IndexType sizesTarget[SCAI_GRID_MAX_DIMENSION];
    IndexType distancesSource[SCAI_GRID_MAX_DIMENSION];
    IndexType distancesTarget[SCAI_GRID_MAX_DIMENSION];

    IndexType dimsSource = other.getDopeVector( offsetSource, sizesSource, distancesSource );
    IndexType dimsTarget = getDopeVector( offsetTarget, sizesTarget, distancesTarget );

    SCAI_ASSERT_EQ_ERROR( dimsSource, dimsTarget, "section dimensions do not match" )

    for ( IndexType i = 0; i < dimsSource; ++i )
    {
        SCAI_ASSERT_EQ_ERROR( sizesSource[i], sizesTarget[i], "size mismatch for section dim = " << i )
    }

    hmemo::ReadAccess<ValueType> rSourceSection( other.mGridVector.getLocalValues() );
    hmemo::WriteAccess<ValueType> wTargetSection( mGridVector.getLocalValues() );

    const ValueType* sourcePtr = rSourceSection.get() + offsetSource;
    ValueType* targetPtr = wTargetSection.get() + offsetTarget;

    if ( dimsSource == 3 )
    {
        for ( IndexType i0 = 0; i0 < sizesSource[0]; ++i0 )
        {
            for ( IndexType i1 = 0; i1 < sizesSource[1]; ++i1 )
            {
                for ( IndexType i2 = 0; i2 < sizesSource[2]; ++i2 )
                {
                    targetPtr[ i0 * distancesTarget[0] + i1 * distancesTarget[1] + i2 * distancesTarget[2]] = 
                        sourcePtr[ i0 * distancesSource[0] + i1 * distancesSource[1] + i2 * distancesSource[2]];
                }
            }
        }
    }
    else
    {
        COMMON_THROWEXCEPTION( "section dims = " << dimsSource << " not supported yet" )
    }

    return *this;
}

template<typename ValueType>
GridSection<ValueType>& GridSection<ValueType>::operator= ( const ValueType other )
{
    GridWriteAccess<ValueType> wTargetSection( mGridVector );

    if ( mNDims ==  3)
    {
        const IndexType* tRange0 = mDimension[0].globalRange;
        const IndexType* tRange1 = mDimension[1].globalRange;
        const IndexType* tRange2 = mDimension[2].globalRange;

        std::cout << "target range 0 = " << tRange0[0] << " - " << tRange0[1] << std::endl;
        std::cout << "target range 1 = " << tRange1[0] << " - " << tRange1[1] << std::endl;
        std::cout << "target range 2 = " << tRange2[0] << " - " << tRange2[1] << std::endl;
    
        for ( IndexType i0 = tRange0[0]; i0 < tRange0[1]; i0++ )
        {
            for ( IndexType i1 = tRange1[0]; i1 < tRange1[1]; i1++ )
            {
                for ( IndexType i2 = tRange2[0]; i2 < tRange2[1]; i2++ )
                {
                    wTargetSection( i0, i1, i2) = other;
                }
            }
        }
    }
    else if ( mNDims == 2 )
    {
        const IndexType* tRange0 = mDimension[0].globalRange;
        const IndexType* tRange1 = mDimension[1].globalRange;

        std::cout << "target range 0 = " << tRange0[0] << " - " << tRange0[1] << std::endl;
        std::cout << "target range 1 = " << tRange1[0] << " - " << tRange1[1] << std::endl;

        for ( IndexType i0 = tRange0[0]; i0 < tRange0[1]; i0++ )
        {
            for ( IndexType i1 = tRange1[0]; i1 < tRange1[1]; i1++ )
            {
                wTargetSection( i0, i1 ) = other;
            }
        }
    }
    else if ( mNDims == 1 )
    {
        const IndexType* tRange0 = mDimension[0].globalRange;

        std::cout << "target range 0 = " << tRange0[0] << " - " << tRange0[1] << std::endl;

        for ( IndexType i0 = tRange0[0]; i0 < tRange0[1]; i0++ )
        {
            wTargetSection( i0 ) = other;
        }
    }
    else
    {
        COMMON_THROWEXCEPTION( "only 2,3 dims are supported" )
    }

    return *this;
}

template<typename ValueType>
GridSection<ValueType>& GridSection<ValueType>::operator*= ( const ValueType other )
{
    IndexType sectionOffset = 0;
    IndexType sectionSizes[SCAI_GRID_MAX_DIMENSION];
    IndexType sectionDistances[SCAI_GRID_MAX_DIMENSION];

    IndexType sectionDims = getDopeVector( sectionOffset, sectionSizes, sectionDistances );

    SCAI_LOG_ERROR( logger, "section dims = " << sectionDims << ", offset = " << sectionOffset );

    for ( IndexType k = 0; k < sectionDims; ++k )
    {
        SCAI_LOG_ERROR( logger, "size[" << k << "] = " << sectionSizes[k] << ", distance = " << sectionDistances[k] )
    }

    SCAI_ASSERT_EQ_ERROR( 3, mNDims, "other than 3 dims not supported yet" );

    GridWriteAccess<ValueType> wTargetSection( mGridVector );

    const IndexType* tRange0 = mDimension[0].globalRange;
    const IndexType* tRange1 = mDimension[1].globalRange;
    const IndexType* tRange2 = mDimension[2].globalRange;

    std::cout << "target range 0 = " << tRange0[0] << " - " << tRange0[1] << std::endl;
    std::cout << "target range 1 = " << tRange1[0] << " - " << tRange1[1] << std::endl;
    std::cout << "target range 2 = " << tRange2[0] << " - " << tRange2[1] << std::endl;

    for ( IndexType i0 = tRange0[0]; i0 < tRange0[1]; i0++ )
    {
        for ( IndexType i1 = tRange1[0]; i1 < tRange1[1]; i1++ )
        {
            for ( IndexType i2 = tRange2[0]; i2 < tRange2[1]; i2++ )
            {
                wTargetSection( i0, i1, i2) *= other;
            }
        }
    }

    return *this;
}

/* ========================================================================= */
/*       Template instantiations                                             */
/* ========================================================================= */

SCAI_COMMON_INST_CLASS( GridSection, SCAI_ARRAY_TYPES_HOST )

} /* end namespace lama */

} /* end namespace scai */
