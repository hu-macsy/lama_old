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

#include <scai/lama/GridReadAccess.hpp>
#include <scai/lama/GridWriteAccess.hpp>

#include <scai/common/macros/instantiate.hpp>
#include <scai/common/SCAITypes.hpp>

namespace scai
{

namespace lama
{

template<typename ValueType>
GridSection<ValueType>& GridSection<ValueType>::operator= ( const GridSection<ValueType>& other )
{
    SCAI_ASSERT_EQ_ERROR( mNDims, other.mNDims, "section mismatch" );

    SCAI_ASSERT_EQ_ERROR( 3, mNDims, "other than 3 dims not supported yet" );

    // okay: section assignment on same array
    // Todo: consider aliasing, i.e. overlapping section

    GridReadAccess<ValueType> rSourceSection( other.mGridVector );
    GridWriteAccess<ValueType> wTargetSection( mGridVector );

    const IndexType* tRange0 = mDimension[0].globalRange;
    const IndexType* tRange1 = mDimension[1].globalRange;
    const IndexType* tRange2 = mDimension[2].globalRange;

    const IndexType* sRange0 = other.mDimension[0].globalRange;
    const IndexType* sRange1 = other.mDimension[1].globalRange;
    const IndexType* sRange2 = other.mDimension[2].globalRange;

    // for ( IndexType i0 = tRange0[0], j0 = sRange0[0]; i0 < tRange0[1]; i0 += tRange0[2], j0 += sRange0[2] )
    // {
        // for ( IndexType i1 = tRange1[0], j1 = sRange1[0]; i1 < tRange1[1]; i1 += tRange1[2], j1 += sRange1[2] )
        // {
            // for ( IndexType i2 = tRange2[0], j2 = sRange2[0]; i2 < tRange2[1]; i2 += tRange2[2], j2 += sRange2[2] )
            // {

    std::cout << "target range 0 = " << tRange0[0] << " - " << tRange0[1] << std::endl;
    std::cout << "target range 1 = " << tRange1[0] << " - " << tRange1[1] << std::endl;
    std::cout << "target range 2 = " << tRange2[0] << " - " << tRange2[1] << std::endl;

    std::cout << "source range 0 = " << sRange0[0] << " - " << sRange0[1] << std::endl;
    std::cout << "source range 1 = " << sRange1[0] << " - " << sRange1[1] << std::endl;
    std::cout << "source range 2 = " << sRange2[0] << " - " << sRange2[1] << std::endl;

    for ( IndexType i0 = tRange0[0], j0 = sRange0[0]; i0 < tRange0[1]; i0++, j0++ )
    {
        for ( IndexType i1 = tRange1[0], j1 = sRange1[0]; i1 < tRange1[1]; i1++, j1++ )
        {
            for ( IndexType i2 = tRange2[0], j2 = sRange2[0]; i2 < tRange2[1]; i2++, j2++ )
            {
                wTargetSection(i0,i1,i2) = rSourceSection( j0, j1, j2 ); 
            }
        }
    }

    return *this;
}

template<typename ValueType>
GridSection<ValueType>& GridSection<ValueType>::operator= ( const ValueType other )
{
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
                wTargetSection( i0, i1, i2) = other;
            }
        }
    }

    return *this;
}

template<typename ValueType>
GridSection<ValueType>& GridSection<ValueType>::operator*= ( const ValueType other )
{
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
