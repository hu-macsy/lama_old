/**
 * @file GPIMemManager.cpp
 *
 * @license
 * Copyright (c) 2009-2016
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
 * @brief Implementation of GPIMemManager routines.
 * @author Thomas Brandes
 * @date 09.05.2014
 */

#include <scai/dmemo/gpi/GPIMemManager.hpp>
#include <scai/dmemo/gpi/GPIUtils.hpp>

#include <scai/common/macros/assert.hpp>
#include <scai/tracing.hpp>
#include <scai/common/SCAITypes.hpp>

// boost
#include <boost/version.hpp>
#include <boost/interprocess/segment_manager.hpp>
#include <boost/interprocess/managed_shared_memory.hpp>
#include <boost/interprocess/managed_heap_memory.hpp>

#include <cstdlib>

using namespace boost::interprocess;

namespace scai
{

namespace dmemo
{
// Help routine to get the Size of the Shared Memory segment

static size_t getSegmentSize()
{
    size_t size = 1024 * 1024;   // MB

    std::cout << "getSegmentSize -> " << size << " * nMB " << std::endl;

    int nMB = 256;   // default 256 MB

    const char* sizeString = getenv( "SCAI_GPI_SEGMENT_SIZE" );

    if ( sizeString  )
    {
        std::cout << "Environment variable SCAI_GPI_SEGMENT_SIZE = " << sizeString << std::endl;
        int n = sscanf( sizeString, "%d", &nMB );
        std::cout << "Read nMB = " << nMB << ", n = " << n << std::endl;
    }
    else
    {
        std::cout << "Environment variable SCAI_GPI_SEGMENT_SIZE not found" << std::endl;
    }

    std::cout << "getSegmentSize -> " << size << " * " << nMB << std::endl;

    size *= nMB;

    std::cout << "getSegmentSize -> " << size << std::endl;

    return size;
}

SCAI_LOG_DEF_LOGGER( GPIMemManager::logger, "GPIMemManager" )

static gaspi_segment_id_t theSegId = 0;

#if BOOST_VERSION >= 105300
typedef ipcdetail::basic_managed_memory_impl
#else
typedef detail::basic_managed_memory_impl
#endif
<char, rbtree_best_fit<mutex_family>, iset_index> base_t;

class ManagedSegmentMemory : public base_t
{
public:

    ManagedSegmentMemory( gaspi_segment_id_t id, std::size_t size )
    {
        mId = id;

        // allocate GASPI segment

        SCAI_GASPI_CALL( gaspi_segment_create ( mId, size, GASPI_GROUP_ALL, GASPI_BLOCK, GASPI_ALLOC_DEFAULT ) )

        gaspi_pointer_t ptr = NULL;

        SCAI_GASPI_CALL( gaspi_segment_ptr ( mId, &ptr ) )

        SCAI_GASPI_CALL( gaspi_barrier( GASPI_GROUP_ALL, GASPI_BLOCK ) )

        base_t::create_impl( ptr, size );

        mPtr = static_cast<char*>( ptr );

        mSize = size;
    }

    ~ManagedSegmentMemory()
    {
    }

    char* get() const
    {
        return mPtr;
    }

    IndexType getOffset( void* ptr )
    {
        char* dataPtr = static_cast<char*>( ptr );

        IndexType offset = dataPtr - mPtr;

        return offset;
    }

    bool getOffset( gaspi_offset_t& offset, const gaspi_pointer_t ptr )
    {
        char* dataPtr = static_cast<char*>( ptr );

        if ( dataPtr < mPtr )
        {
            offset = 0;
            return false;
        }

        if ( dataPtr >= mPtr + mSize )
        {
            offset = 0;
            return false;
        }

        offset = dataPtr - mPtr;

        return true;
    }

private:

    gaspi_segment_id_t mId;
    char* mPtr;
    size_t mSize;
};

static ManagedSegmentMemory* theSegmentManager = NULL;

void GPIMemManager::getSegmentData( gaspi_segment_id_t& id, gaspi_pointer_t& ptr, gaspi_offset_t& offset, const int size )
{
    SCAI_REGION( "GPIMemManager.getSegmentData" )

    if ( theSegmentManager == NULL )
    {
        size_t segmentSize = getSegmentSize();

        SCAI_LOG_WARN( logger, "Create GASPI segment of size " << segmentSize )

        theSegmentManager = new ManagedSegmentMemory( theSegId, segmentSize );
    }

    id = theSegId;

    ptr = 0;

    try
    {
        ptr = theSegmentManager->allocate( size );
    }
    catch ( std::exception& ex )
    {
        SCAI_LOG_ERROR( logger, "allocate( " << size << " ) failed, " << ex.what() )
    }

    if ( ptr == 0 )
    {
        COMMON_THROWEXCEPTION( "Serious allocation error, size = " << size )
    }

    offset = theSegmentManager->getOffset( ptr );

    SCAI_LOG_DEBUG( logger, "getSegmentData( size = " << size << " ) =>, offset = " << offset << ", ptr = " << ptr )
}

bool GPIMemManager::findSegment( gaspi_segment_id_t& id, gaspi_offset_t& offset, const gaspi_pointer_t ptr )
{
    if ( theSegmentManager == NULL )
    {
        return false;
    }

    id = theSegId;

    return theSegmentManager->getOffset( offset, ptr );
}

void GPIMemManager::releaseSegmentData( const gaspi_segment_id_t id, const gaspi_offset_t offset )
{
    SCAI_ASSERT_ERROR( theSegmentManager != NULL, "release segment data, but no manager" )

    char* dataPtr = theSegmentManager->get() + offset;

    SCAI_LOG_DEBUG( logger, "releaseSegmentData, id = " << static_cast<int>( id ) << ", offset = " << offset )

    theSegmentManager->deallocate( dataPtr );
}

void GPIMemManager::freeAll( )
{
    if ( theSegmentManager == NULL )
    {
        return;
    }

    delete theSegmentManager;

    theSegmentManager = NULL;
}

} // namespace dmemo

} // namespace scai

