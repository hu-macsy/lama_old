/**
 * @file GeneralDistribution.cpp
 *
 * @license
 * Copyright (c) 2011
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 * @endlicense
 *
 * @brief GeneralDistribution.cpp
 * @author brandes
 * @date 25.02.2011
 * $Id$
 */

// hpp
#include <lama/distribution/GeneralDistribution.hpp>

// assert
#include <lama/exception/LAMAAssert.hpp>

#include <algorithm>
#include <functional>

//#include <fstream>

#define MASTER 0

namespace lama
{

LAMA_LOG_DEF_LOGGER( GeneralDistribution::logger, "Distribution.General" );

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
        LAMA_ASSERT( 0 <= *it && *it < mGlobalSize,
                     *it << " is illegal index for general distribution of size " << mGlobalSize );
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
    LAMA_LOG_INFO( logger, "~GeneralDistribution" );
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
        LAMA_ASSERT( displ[ parts ] == mGlobalSize, "sum of local rows is not global size" );
    }

    // gather global indices of local rows
    std::vector<IndexType> rows( mGlobalSize );

    mCommunicator->gather( &rows[0], numMyRows, MASTER, &mLocal2Global[0], &numRows[0] );

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

void GeneralDistribution::printDistributionVector( std::string /*name*/) const
{
//    IndexType myRank = mCommunicator->getRank();
    IndexType parts = mCommunicator->getSize();

    // gather number of local rows
    IndexType numMyRows = static_cast<IndexType>( mLocal2Global.size() );
    std::vector<IndexType> numRows( parts );
    mCommunicator->gather( &numRows[0], 1, MASTER, &numMyRows );

    // gather global indices of local rows
    std::vector<IndexType> rows( mGlobalSize );

    mCommunicator->gather( &rows[0], numMyRows, MASTER, &mLocal2Global[0], &numRows[0] );

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

}
