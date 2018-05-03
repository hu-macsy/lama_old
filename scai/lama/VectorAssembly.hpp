/**
 * @file VectorAssembly.hpp
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

#include <scai/dmemo/Distribution.hpp>
#include <scai/hmemo/HArray.hpp>
#include <scai/common/Printable.hpp>

#include <scai/logging.hpp>

#include <vector>

namespace scai
{

namespace lama
{

/** This classs allows to assembly vector entries by different processors. Each processor can 
 *  add vector entries by global index and the value. 
 *
 *  An object stands for a set of collected vector entries and can be converted into a vector 
 *  with a given distribution.
 */
template<typename ValueType>
class VectorAssembly : public common::Printable
{

public:

    /** Construct an empty vector assembly
     *
     *  @param comm specifies the processor set onto which the data will be assembled
     *
     *  By default, elements are assembled by all processors of the current communicator. 
     *  The communicator might be set explicitly to NoCommunicator to indicate that elements
     *  are assembled in a repicated way.
     */
    VectorAssembly( const dmemo::CommunicatorPtr comm = dmemo::Communicator::getCommunicatorPtr() );

    /** Destructor */

    ~VectorAssembly();

    /** Might be used for optimization to indicate how many elements might be added by this processor. */

    void reserve( const IndexType n );

    /** Add a vector element with global coordinates, called by a single processor */

    void push( const IndexType i, const ValueType val );

    /**
     *  Override Printable::writeAt 
     */
    virtual void writeAt( std::ostream& stream ) const;

    /**
     *  Query the minimal global size of the vector that might set/add this assembly.
     */
    IndexType getSize() const;

    /** Query the number of entries collected in the assembly so far.  */

    IndexType getNumValues() const;

    /**
     *  @brief Get the assembled sparse data localized for a given distribution
     *
     *  @param[out] ia     local indexes of assembled data owned by this processor correspoding dist
     *  @param[out] values the corresponding values at the positions given by indexes
     *  @param[in]  dist   specifies the ownership how to redistribute the assembled data.
     *  
     *  Be careful: the entries are not sorted and might contain multiple values.
     *
     */
    void buildLocalData( 
        hmemo::HArray<IndexType>& ia, 
        hmemo::HArray<ValueType>& values, 
        const dmemo::Distribution& dist ) const;

    /**
     *  @brief Get the assembled sparse data replicated for all processors
     */
    void buildGlobalData( 
        hmemo::HArray<IndexType>& ia, 
        hmemo::HArray<ValueType>& values,
        const IndexType n ) const;

    /**
     *  @brief Query the communicator (processor set) of assembled data
     */
    const dmemo::Communicator& getCommunicator() const;

private:

    SCAI_LOG_DECL_STATIC_LOGGER( logger )

    /** Resort sparse vector data according to the ownership of the indexes */

    void exchangeCOO(                      
        hmemo::HArray<IndexType>& outIA,
        hmemo::HArray<ValueType>& outValues,
        const hmemo::HArray<IndexType> inIA,
        const hmemo::HArray<ValueType> inValues,
        const dmemo::Distribution& dist ) const;

    /** Check for correct indexes of the assembled entries */

    void checkLegalIndexes( const IndexType size ) const;

    // for pushing globally assembled data we use the C++ vector class

    std::vector<IndexType> mIA;
    std::vector<ValueType> mValues;

    dmemo::CommunicatorPtr mComm;
};

/* ================================================================================ */
/*   Implementation of inline methods                                               */
/* ================================================================================ */

template<typename ValueType>
void VectorAssembly<ValueType>::push( const IndexType i, const ValueType val )
{
    mIA.push_back( i );
    mValues.push_back( val );
}

template<typename ValueType>
const dmemo::Communicator& VectorAssembly<ValueType>::getCommunicator() const
{
    return *mComm;
}

}

}
