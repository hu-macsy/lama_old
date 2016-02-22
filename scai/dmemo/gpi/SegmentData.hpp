/**
 * @file SegmentData.hpp
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
 * @brief Class to deal with dynamic arrays in GASPI segments
 * @author Thomas Brandes
 * @date 09.05.2014
 * @since 1.1.0
 */

#pragma once

#include <scai/hmemo/Access.hpp>

#include <scai/common/SCAITypes.hpp>
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

template<typename T>
class SegmentData : public hmemo::Access
{
public:

    /** Constructor for an array of type T with size elements. */

    SegmentData( const GPICommunicator* comm, const IndexType size );

    /** Constructor for an array of type T with available data */

    SegmentData( const GPICommunicator* comm, const IndexType size, T* data );

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
    void assign( const T values[], const IndexType n );

    /** Copy data from SegmentData array to other memory
     *
     *  @param[out] values points to data to be filled
     *  @param[in] n  number of data values to be copied
     */
    void copyTo( T values[], const IndexType n ) const;

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
        return mOffsetBytes / sizeof( T );
    }

    // return typed pointer to the segment data

    T* get() 
    {
        return mData;
    }

    /** Get const typed pointer to the segment data. */

    const T* get() const
    {
        return mData;
    }

    /** Operator [] allows access to single values of the array. */

    T& operator[]( const IndexType i )
    {
        return mData[i];
    }

    /** Useful info about object written into a stream. */

    virtual void writeAt( std::ostream& stream ) const;

private:

    void reserve( const IndexType size );

    const GPICommunicator* mComm;

    gaspi_pointer_t    mPtr;    // pointer to the allocated segment data
    T*                 mData;   // might be aligned, i.e. mPtr + <align_bytes>

    gaspi_segment_id_t mId;
    gaspi_offset_t     mOffsetBytes;   // mPtr = segmentPtr<mId> + mOffsetBytes

    IndexType      mSize;    // might be used for checks, infos

    bool mReleaseFlag;    // if true segment data must be released

    SCAI_LOG_DECL_STATIC_LOGGER( logger )
};

} // namespace dmemo

} // namespace scai
