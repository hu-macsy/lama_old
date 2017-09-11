/**
 * @file MatrixAssemblyAccess.hpp
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
 * @brief Access to a matrix to add matrix elements
 * @author Thomas Brandes
 * @date 07.09.2017
 */

#pragma once 

#include <scai/lama.hpp>

#include <scai/lama/matrix/Matrix.hpp>

#include <vector>

namespace scai
{

namespace lama
{

/** This classs allows to assembly matrix entries by different processors. Each processor can 
 *  add matrix entries by global coordinates. Non-local entries can also be pushed as these
 *  elements will be redistributed when the access is released.
 */

template<typename ValueType>
class MatrixAssemblyAccess 
{

public:

    /** Construct an access */

    MatrixAssemblyAccess( Matrix& matrix, const common::binary::BinaryOp op = common::binary::COPY );

    /** Destructor of the access, also inserts the assembled entries into the matrix. */

    ~MatrixAssemblyAccess()
    {
        if ( !mIsReleased )
        {
            release();
        }
    }

    /** Might be used for optimization to indicate how many elements might be added by this processor. */

    void reserve( const IndexType n )
    {
        mIA.reserve( n );
        mJA.reserve( n );
        mValues.reserve( n );
    }

    /** Add a matrix element with global coordinates */

    void push( const IndexType i, const IndexType j, const ValueType val )
    {
        SCAI_ASSERT_VALID_INDEX_DEBUG( i, mMatrix.getNumRows(), "illegal row index pushed" );
        SCAI_ASSERT_VALID_INDEX_DEBUG( j, mMatrix.getNumColumns(), "illegal column index pushed" );

        mIA.push_back( i );
        mJA.push_back( j );
        mValues.push_back( val );

        SCAI_LOG_TRACE( logger, mMatrix.getRowDistribution().getCommunicator() << ": pushed " << val << " @ ( " << i << ", " << j << " )" )
    }

    /** Release the assembly access, all pushed entries are now transferred to owning processors and added. */

    void release();

private:

    SCAI_LOG_DECL_STATIC_LOGGER( logger )

    /** Translate global indexes to local indexes */

    void global2local( hmemo::HArray<IndexType>& ia,
                       const dmemo::Distribution& dist );

    /** Resort COO data according to the ownership of the row indexes */

    void exchangeCOO(                      
        hmemo::HArray<IndexType>& outIA,
        hmemo::HArray<IndexType>& outJA,
        hmemo::HArray<ValueType>& outValues,
        const hmemo::HArray<IndexType> inIA,
        const hmemo::HArray<IndexType> inJA,
        const hmemo::HArray<ValueType> inValues,
        const dmemo::Distribution& dist );

    Matrix& mMatrix;

    // for pushing the assembled data we use the C++ vector class

    std::vector<IndexType> mIA;
    std::vector<IndexType> mJA;
    std::vector<ValueType> mValues;

    bool mIsReleased;

    common::binary::BinaryOp mOp;   // specifies how to combine with available entries
};

}

}
