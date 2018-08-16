/**
 * @file JoinedDistribution.hpp
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
 * @brief Distribution class that stands for the concatenation of two distributions
 * @author Thomas Brandes
 * @date 27.07.2017
 */

#pragma once

#include <scai/dmemo/Distribution.hpp>

namespace scai

{

namespace dmemo

{

/** An object of this class stands for a distribution of joined vectors
 *  or of joined matrices (first dimension)
 */
class JoinedDistribution : public Distribution
{

public:

    /** Constructor concatenates two distributions 
     * 
     *  @param[in] d1 distribution of first part
     *  @param[in] d2 distribution of second part
     * 
     *  As distributions are always shared pointers, the use is safe regarding 
     *  the lifetime of the input arguments.
     */
    JoinedDistribution ( DistributionPtr d1, DistributionPtr d2 ) :

        Distribution( d1->getGlobalSize() + d2->getGlobalSize() ),
        mD1( d1 ),
        mD2( d2 )

    {
        // no restrictions on the distribution, all can be joined
    }

    ~JoinedDistribution()
    {
    }

 	virtual const char* getKind() const
    {
        return "JOINED";
    }

  	virtual bool isLocal( IndexType globalIndex ) const
    {
        if ( globalIndex < mD1->getGlobalSize() )
        {
            return mD1->isLocal( globalIndex );
        }
        else
        {
            return mD2->isLocal( globalIndex - mD1->getGlobalSize() );
        }
    }

  	virtual IndexType getLocalSize() const
    {
        return mD1->getLocalSize() + mD2->getLocalSize();
    }

  	virtual IndexType local2global( IndexType localIndex ) const
    {
        const IndexType localSize1 = mD1->getLocalSize();

        if ( localIndex < localSize1 )
        {
            // look up in first distribution
            return mD1->local2global( localSize1 );
        }
        else
        {
            // look up in second distribution
            return mD1->getGlobalSize() + mD2->local2global( localIndex - localSize1 );
        }
    }

  	virtual IndexType global2local( IndexType globalIndex ) const
    {
        const IndexType globalSize1 = mD1->getGlobalSize();
  
        if ( globalIndex < globalSize1 )
        {
            return global2local( globalIndex );
        }
        else
        {
            IndexType localIndex = mD2->global2local( globalIndex - globalSize1 );

            if ( localIndex != invalidIndex )
            {
                localIndex += mD1->getLocalSize();
            }

            return localIndex;
        }
    }

  	virtual IndexType getBlockDistributionSize() const
    {
        IndexType blockSize1 = mD1->getBlockDistributionSize();
        IndexType blockSize2 = mD2->getBlockDistributionSize();

        if ( blockSize1 != invalidIndex && blockSize2 != invalidIndex )
        {
            return blockSize1 + blockSize2;
        }
        else
        {
            return invalidIndex;
        }
    }

  	virtual void enableAnyAddressing() const
    {
        mD1->enableAnyAddressing();
        mD2->enableAnyAddressing();
    }

    virtual bool hasAnyAddressing() const
    {
        return mD1->hasAnyAddressing() && mD2->hasAnyAddressing();
    }

  	virtual IndexType getAnyLocalSize( PartitionId p ) const
    {
        return mD1->getAnyLocalSize( p ) + mD2->getAnyLocalSize( p );
    }

  	virtual IndexType getAnyOwner(IndexType) const
    {
        COMMON_THROWEXCEPTION( "getAnyOwner not availalbe yet" )
        return 0;
    }

  	virtual IndexType getAnyLocalIndex(IndexType , PartitionId) const
    {
        COMMON_THROWEXCEPTION( "getAnyLocalIndex not availalbe yet" )
        return 0;
    }

  	virtual IndexType getAnyGlobalIndex(IndexType, PartitionId) const
    {
        COMMON_THROWEXCEPTION( "getAnyGlobalIndex not availalbe yet" )
        return 0;
    }

  	virtual bool isEqual( const Distribution& other ) const
    {
        if ( strcmp( other.getKind(), getKind() ) != 0 )
        {
            return false;
        }

        const JoinedDistribution& jOther = static_cast<const JoinedDistribution&>( other );
       
        return mD1->isEqual( *jOther.mD1 ) && mD2->isEqual( *jOther.mD2 );
    }

private:

    DistributionPtr mD1;
    DistributionPtr mD2;
};

}

}
