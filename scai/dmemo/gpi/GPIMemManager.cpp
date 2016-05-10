/**
 * @file GPIMemManager.cpp
 *
 * @license
 * Copyright (c) 2009-2013
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
 * @brief Implementation of GPIMemManager routines.
 * @author Thomas Brandes
 * @date 09.05.2014
 * @since 1.1.0
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

