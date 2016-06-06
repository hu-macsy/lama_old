/**
 * @file GeneralDistribution.cpp
 *
 * @license
 * Copyright (c) 2009-2016
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the Library of Accelerated Math Applications (LAMA).
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
 * @endlicense
 *
 * @brief GeneralDistribution.cpp
 * @author brandes
 * @date 25.02.2011
 */

// hpp
#include <scai/dmemo/GeneralDistribution.hpp>

// internal scai libraries
#include <scai/common/macros/assert.hpp>

// std
#include <algorithm>
#include <functional>

#define MASTER 0

namespace scai
{

namespace dmemo
{

SCAI_LOG_DEF_LOGGER( GeneralDistribution::logger, "Distribution.General" )

const char GeneralDistribution::theCreateValue[] = "GENERAL";

GeneralDistribution::GeneralDistribution(
    const IndexType globalSize,
    const std::vector<IndexType>& localIndexes,
    const CommunicatorPtr communicator )
    : Distribution( globalSize, communicator ), mLocal2Global( localIndexes )
{
    std::vector<IndexType>::const_iterator end = mLocal2Global.end();
    std::vector<IndexType>::const_iterator begin = mLocal2Global.begin();

    for ( std::vector<IndexType>::const_iterator it = begin; it != end; ++it )
    {
        IndexType i = static_cast<IndexType>( std::distance( begin, it ) );
        SCAI_ASSERT( 0 <= *it && *it < mGlobalSize,
                     *it << " is illegal index for general distribution of size " << mGlobalSize )
        mGlobal2Local[ *it] = i;
    }
}

GeneralDistribution::GeneralDistribution(
    const std::vector<IndexType>& row2Partition,
    const IndexType globalSize,
    const CommunicatorPtr communicator )
    : Distribution( globalSize, communicator )
{
    IndexType myRank = mCommunicator->getRank();
    IndexType parts = mCommunicator->getSize();
    IndexType partSize = row2Partition.size();
    IndexType numMyRows = 0;

    // // gather number of local rows
    // IndexType numMyRows = static_cast<IndexType>( mLocal2Global.size() );
    // SCAI_ASSERT(mGlobalSize == partSize, "partition size " << partSize << " is not equal to global size " << mGlobalSize);

    std::vector<IndexType> displ;
    std::vector<IndexType> curpos;
    std::vector<IndexType> rows;

    if ( myRank == MASTER )
    {
        SCAI_ASSERT( mGlobalSize == partSize,
                     "partition size " << partSize << " is not equal to global size " << mGlobalSize )
        displ.resize( parts + 1 );

        for ( IndexType i = 0; i < mGlobalSize; ++i )
        {
            SCAI_ASSERT( row2Partition[i] < parts, "invalid partition id at position" << i )
            displ[row2Partition[i] + 1]++;
        }
    }
    else
    {
        displ.resize( 2 );
    }

    // scatter partition sizes
    mCommunicator->scatter( &numMyRows, 1, MASTER, &displ[1] );

    if ( myRank == MASTER )
    {
        rows.resize( mGlobalSize );
        curpos.resize( parts );

        for ( IndexType i = 1; i < parts; i++ )
        {
            displ[i + 1] += displ[i];
            curpos[i] = 0;
        }

        SCAI_ASSERT( displ[parts] == mGlobalSize, "sum of local rows is not global size" )

        for ( IndexType i = 0; i < mGlobalSize; ++i )
        {
            IndexType partition = row2Partition[i];
            IndexType position = displ[partition] + curpos[partition]++;
            rows[position] = i;
        }

        for ( IndexType i = 0; i < parts; ++i )
        {
            SCAI_ASSERT( displ[i] + curpos[i] == displ[i + 1],
                         "partition " << i << "  size mismatch, expected " << displ[i + 1] - displ[i] << " actual " << curpos[i] )
        }
    }
    else
    {
        rows.resize( 1 );
        curpos.resize( 1 );
    }

    // scatter global indices of local rows
    mLocal2Global.resize( numMyRows );
    mCommunicator->scatterV( &mLocal2Global[0], numMyRows, MASTER, &rows[0], &curpos[0] );

    // Compute Global2Local
    std::vector<IndexType>::const_iterator end = mLocal2Global.end();
    std::vector<IndexType>::const_iterator begin = mLocal2Global.begin();

    for ( std::vector<IndexType>::const_iterator it = begin; it != end; ++it )
    {
        IndexType i = static_cast<IndexType>( std::distance( begin, it ) );
        SCAI_ASSERT( 0 <= *it && *it < mGlobalSize,
                     *it << " is illegal index for general distribution of size " << mGlobalSize )
        mGlobal2Local[ *it] = i;
    }
}

GeneralDistribution::GeneralDistribution( const Distribution& other )
    : Distribution( other.getGlobalSize(), other.getCommunicatorPtr() ), mLocal2Global(
          other.getLocalSize() )
{
    for ( IndexType i = 0; i < getGlobalSize(); ++i )
    {
        if ( other.isLocal( i ) )
        {
            IndexType localIndex = other.global2local( i );
            mGlobal2Local[i] = localIndex;
            mLocal2Global[localIndex] = i;
        }
    }
}

GeneralDistribution::GeneralDistribution( const IndexType globalSize, const CommunicatorPtr communicator )
    : Distribution( globalSize, communicator )
{
}

GeneralDistribution::~GeneralDistribution()
{
    SCAI_LOG_INFO( logger, "~GeneralDistribution" )
}

bool GeneralDistribution::isLocal( const IndexType index ) const
{
    return mGlobal2Local.find( index ) != mGlobal2Local.end();
}

IndexType GeneralDistribution::getLocalSize() const
{
    return static_cast<IndexType>( mLocal2Global.size() );
}

std::vector<IndexType>& GeneralDistribution::getLocalRows()
{
    return mLocal2Global;
}

IndexType GeneralDistribution::local2global( const IndexType localIndex ) const
{
    return mLocal2Global[localIndex];
}

IndexType GeneralDistribution::global2local( const IndexType globalIndex ) const
{
    const Global2LocalMapType::const_iterator elem = mGlobal2Local.find( globalIndex );

    if ( elem == mGlobal2Local.end() )
    {
        return nIndex;
    }

    return elem->second;
}

bool GeneralDistribution::isEqual( const Distribution& other ) const
{
    return this == &other;
}

void GeneralDistribution::writeAt( std::ostream& stream ) const
{
    // write identification of this object

    stream << "GeneralDistribution( size = " << mLocal2Global.size() << " of " << mGlobalSize << ", comm = "
           << *mCommunicator << " )";
}

void GeneralDistribution::getDistributionVector( std::vector<IndexType>& row2Partition ) const
{
    IndexType myRank = mCommunicator->getRank();
    IndexType parts = mCommunicator->getSize();

    // gather number of local rows
    IndexType numMyRows = static_cast<IndexType>( mLocal2Global.size() );
    std::vector<IndexType> numRows( parts );
    mCommunicator->gather( &numRows[0], 1, MASTER, &numMyRows );

    std::vector<IndexType> displ;

    if ( myRank == MASTER )
    {
        displ.reserve( parts + 1 );

        IndexType displacement = 0;

        for ( IndexType i = 0; i < parts; i++ )
        {
            displ[i] = displacement;
            displacement += numRows[i];
        }

        displ[parts] = displacement;
        SCAI_ASSERT( displ[parts] == mGlobalSize, "sum of local rows is not global size" )
    }

    // gather global indices of local rows
    std::vector<IndexType> rows( mGlobalSize );

    mCommunicator->gatherV( &rows[0], numMyRows, MASTER, &mLocal2Global[0], &numRows[0] );

    // build mapping row 2 partition
    if ( myRank == MASTER )
    {

        // for testing: init
        for ( IndexType i = 0; i < mGlobalSize; ++i )
        {
            row2Partition[i] = -1;
        }

        for ( IndexType i = 0; i < parts; ++i )
        {
            for ( IndexType j = displ[i]; j < displ[i + 1]; ++j )
            {
                row2Partition[rows[j]] = i;
            }
        }
    }
}

void GeneralDistribution::printDistributionVector( std::string /*name*/ ) const
{
//    IndexType myRank = mCommunicator->getRank();
    IndexType parts = mCommunicator->getSize();

    // gather number of local rows
    IndexType numMyRows = static_cast<IndexType>( mLocal2Global.size() );
    std::vector<IndexType> numRows( parts );
    mCommunicator->gather( &numRows[0], 1, MASTER, &numMyRows );

    // gather global indices of local rows
    std::vector<IndexType> rows( mGlobalSize );

    mCommunicator->gatherV( &rows[0], numMyRows, MASTER, &mLocal2Global[0], &numRows[0] );

    std::vector<IndexType> row2Partition( mGlobalSize );
    getDistributionVector( row2Partition );
    // build mapping row 2 partition
//    if(myRank == MASTER)
//    {
//
//        std::ofstream file;
//        file.open((name + ".part").c_str());
    // print row - partition mapping
//        std::cout << "partitionVector ";
//        for(IndexType i = 0; i < mGlobalSize; ++i)
//        {
//            file << row2Partition[ i ] << std::endl;
//            std::cout << row2Partition[ i ] << " ";
//        }
//        std::cout << std::endl;
//        file.close();
//    }
}

} /* end namespace dmemo */

} /* end namespace scai */
