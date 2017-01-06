/**
 * @file SegmentData.hpp
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
 * @brief Class to deal with dynamic arrays in GASPI segments
 * @author Thomas Brandes
 * @date 09.05.2014
 */

#pragma once

#include <scai/hmemo/Access.hpp>

#include <scai/common/SCAITypes.hpp>
#include <scai/common/ScalarType.hpp>
#include <scai/dmemo/gpi/GPIMemManager.hpp>

// logging
#include <scai/logging.hpp>

namespace scai
{

namespace dmemo
{

class GPICommunicator;

/** Structure that behaves like a scoped array that can be used for communication.
 *  The array is released/freed by its destructor.
 *
 *  It is derived from Access so it can be used in SyncToken.
 *
 *  @tparam T specifies the data type for the scoped array
 *
 *  Note: In previous versions it was necessary that all processors allocate
 *        SegmentData with the same sizes and in the same order. Otherwise
 *        it could not be sure that the data comes from the same segment id.
 *
 *        This restriction is no more necessary as we have now one big segment
 *        that is managed here itself.
 */

class SegmentData : public hmemo::Access
{
public:

    /** Constructor for an array of type T with size elements. */

    SegmentData( const common::scalar::ScalarType stype, const GPICommunicator* comm, const IndexType size );

    /** Constructor for an array of type T with available data */

    SegmentData( const common::scalar::ScalarType stype, const GPICommunicator* comm, const IndexType size, void* data );

    /** Destructor will release the used memory for other usage. */

    ~SegmentData()
    {
        release();
    }

    /** This routine makes the struct to a BaseAccess.
     *
     *  Segment will be freed for next usage.
     */
    virtual void release();

    /** Assign data to the SegmentData array.
     *
     *  @param[in] values  points the data to be assigned
     *  @param[in] n  number of data values to be assigned
     *
     *  Note: assert n <= mSize
     */
    void assign( const void* values, const IndexType n );

    /** Copy data from SegmentData array to other memory
     *
     *  @param[out] values points to data to be filled
     *  @param[in] n  number of data values to be copied
     */
    void copyTo( void* values, const IndexType n ) const;

    /** Query the segment id to use routines for remote write/read. */

    gaspi_segment_id_t getID() const
    {
        return mId;
    }

    /** Query the offset in bytes to use routines for remote write/read. */

    gaspi_offset_t getOffsetBytes() const
    {
        return mOffsetBytes;
    }

    IndexType getOffset() const
    {
        // get rid of this routine as it does not work for unaligned data

        return mOffsetBytes / mTypeSize;
    }

    IndexType typeSize() const
    {
        return mTypeSize;
    }

    common::scalar::ScalarType scalarType() const
    {
        return mScalarType;
    }

    // return typed pointer to the segment data

    void* get()
    {
        return mData;
    }

    void* get( IndexType offset )
    {
        char* ptr = reinterpret_cast<char*>( mData );
        return ptr + offset * mTypeSize;
    }

    /** Get const typed pointer to the segment data. */

    const void* get() const
    {
        return mData;
    }

    const void* get( IndexType offset ) const
    {
        const char* ptr = reinterpret_cast<const char*>( mData );
        return ptr + offset * mTypeSize;
    }

    /** Operator [] allows access to single values of the array. */

    /*
    T& operator[]( const IndexType i )
    {
        return mData[i];
    }
    */

    /** Useful info about object written into a stream. */

    virtual void writeAt( std::ostream& stream ) const;

private:

    void reserve( const IndexType size );

    const GPICommunicator* mComm;

    common::scalar::ScalarType mScalarType;
 
    IndexType mTypeSize;

    gaspi_pointer_t    mPtr;    // pointer to the allocated segment data
    void*              mData;   // might be aligned, i.e. mPtr + <align_bytes>

    gaspi_segment_id_t mId;
    gaspi_offset_t     mOffsetBytes;   // mPtr = segmentPtr<mId> + mOffsetBytes

    IndexType      mSize;    // might be used for checks, infos

    bool mReleaseFlag;    // if true segment data must be released

    SCAI_LOG_DECL_STATIC_LOGGER( logger )
};

} // namespace dmemo

} // namespace scai
