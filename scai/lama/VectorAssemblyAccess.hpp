/**
 * @file VectorAssemblyAccess.hpp
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
 * @brief Distributed access to a vector to set elements individually
 * @author Thomas Brandes
 * @date 26.09.2017
 */
#pragma once 

#include <scai/lama.hpp>

#include <scai/lama/_Vector.hpp>

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
class VectorAssemblyAccess 
{

public:

    /** Construct an assembly access to a given vector
     *
     *  @param[in] vector is the vector for which data is assembled
     *  @param[in] op specifies how to deal with entries at same positions
     */

    VectorAssemblyAccess( _Vector& vector, const common::BinaryOp op = common::BinaryOp::COPY );

    /** Destructor of the access, also communicates and inserts the assembled entries into the vector. */

    ~VectorAssemblyAccess()
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
        mValues.reserve( n );
    }

    /** Add a vector element with global coordinates, called by a single processor */

    void push( const IndexType i, const ValueType val )
    {
        SCAI_ASSERT_VALID_INDEX_DEBUG( i, mVector.size(), "illegal index pushed" );

        mIA.push_back( i );
        mValues.push_back( val );

        SCAI_LOG_TRACE( logger, mVector.getDistribution().getCommunicator() << ": pushed " << val << " @ ( " << i << " )" )
    }

    /** Add a vector element with global coordinates, called by all processors 
     *  (only the owner process(es) will push it)
     */

    void pushReplicated( const IndexType i, const ValueType val )
    {
        SCAI_ASSERT_VALID_INDEX_DEBUG( i, mVector.size(), "illegal index pushed" );

        const IndexType localI = mVector.getDistribution().global2local( i );

        if ( localI != nIndex )
        {
            mLocalIA.push_back( localI );
            mLocalValues.push_back( val );
        
            SCAI_LOG_TRACE( logger, mVector.getDistribution().getCommunicator() << ": pushed " << val << " @ ( " << i << " )" )
        }
    }

    /** Release the assembly access, all pushed entries are now transferred to owning processors and added. */

    void release();

private:

    SCAI_LOG_DECL_STATIC_LOGGER( logger )

    /** Resort sparse vector data according to the ownership of the indexes */

    void exchangeCOO(                      
        hmemo::HArray<IndexType>& outIA,
        hmemo::HArray<ValueType>& outValues,
        const hmemo::HArray<IndexType> inIA,
        const hmemo::HArray<ValueType> inValues,
        const dmemo::Distribution& dist );

    /** Shifts assembled data so that every processor sees it */

    void shiftAssembledData( 
        const hmemo::HArray<IndexType>& myIA, 
        const hmemo::HArray<ValueType>& myValues );

    _Vector& mVector;

    // for pushing globally assembled data we use the C++ vector class

    std::vector<IndexType> mIA;
    std::vector<ValueType> mValues;

    // locally assembled data in own structures, no more communication needed

    std::vector<IndexType> mLocalIA;
    std::vector<ValueType> mLocalValues;

    bool mIsReleased;

    common::BinaryOp mOp;   // specifies how to combine with available entries
};

}

}
