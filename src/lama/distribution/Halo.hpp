/**
 * @file Halo.hpp
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
 * @brief Halo.hpp
 * @author Thomas Brandes
 * @date 23.02.2011
 * @since 1.0.0
 */
#ifndef LAMA_HALO_HPP_
#define LAMA_HALO_HPP_

// for dll_import
#include <lama/config.hpp>

// base classes
#include <lama/Printable.hpp>

// others
#include <lama/LAMAArray.hpp>
#include <lama/CommunicationPlan.hpp>

#include <lama/exception/LAMAAssert.hpp>

#include <map>

namespace lama
{

/** The halo is an internal data structure that describes the
 *  exchange of non-local values completely.
 *
 *  It is build by the required (global) indexes to set up
 *  communication plans to receive required data and to send
 *  data provided for other partitions.
 */

class LAMA_DLL_IMPORTEXPORT Halo: public Printable
{
    friend class HaloBuilder;

public:

    /** Constructor of a new 'empty' halo. */

    Halo();

    /** Copy constructor. */

    Halo( const Halo& halo );

    virtual ~Halo();

    /** Clear the halo for zero matrix. */

    void clear();

    Halo& operator=( const Halo& other );

    inline const CommunicationPlan& getRequiredPlan() const;

    inline const CommunicationPlan& getProvidesPlan() const;

    inline IndexType global2halo( const IndexType globalIndex ) const;

    inline const LAMAArray<IndexType>& getProvidesIndexes() const;

    inline const LAMAArray<IndexType>& getRequiredIndexes() const;

    /** Query the size for a halo to be allocated */

    inline IndexType getHaloSize() const;

    /** If a halo is empty, no communication is needed for this partition.
     Be careful: getHaloSize() == 0 implies that no indexes are required
     but it might be possible that this partition has to provide values
     */

    inline bool isEmpty() const;

    inline const std::map<IndexType,IndexType>& getMap() const
    {
        return mGlobal2Halo;
    }

    virtual void writeAt( std::ostream& stream ) const;

protected:

    inline void setGlobal2Halo( IndexType globalIndex, IndexType haloIndex );

private:

    CommunicationPlan mRequiredPlan;
    CommunicationPlan mProvidesPlan;

    // Indexes for required values and values to provide are stored in LAMAArrays
    // so they might be used in different contexts, especially also on GPU

    LAMAArray<IndexType> mRequiredIndexes;
    LAMAArray<IndexType> mProvidesIndexes;

    std::map<IndexType,IndexType> mGlobal2Halo;

    LAMA_LOG_DECL_STATIC_LOGGER( logger )
};

const CommunicationPlan& Halo::getRequiredPlan() const
{
    return mRequiredPlan;
}

const CommunicationPlan& Halo::getProvidesPlan() const
{
    return mProvidesPlan;
}

const LAMAArray<IndexType>& Halo::getProvidesIndexes() const
{
    return mProvidesIndexes;
}

const LAMAArray<IndexType>& Halo::getRequiredIndexes() const
{
    return mRequiredIndexes;
}

void Halo::setGlobal2Halo( IndexType globalIndex, IndexType haloIndex )
{
    LAMA_ASSERT_DEBUG( 0 <= haloIndex && haloIndex < getHaloSize(),
                       "illegal halo index " << haloIndex << ", halo size = " << getHaloSize() )
    mGlobal2Halo[globalIndex] = haloIndex;
}

IndexType Halo::global2halo( const IndexType globalIndex ) const
{
    const std::map<IndexType,IndexType>::const_iterator elem = mGlobal2Halo.find( globalIndex );

    if ( elem == mGlobal2Halo.end() )
    {
        return nIndex;
    }

    return elem->second;
}

IndexType Halo::getHaloSize() const
{
    return mRequiredPlan.totalQuantity();
}

bool Halo::isEmpty() const
{
    return ( mRequiredPlan.totalQuantity() == 0 ) && ( mProvidesPlan.totalQuantity() == 0 );
}

}

#endif // LAMA_HALO_HPP_
