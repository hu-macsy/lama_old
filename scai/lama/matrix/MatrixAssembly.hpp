/**
 * @file MatrixAssembly.hpp
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
 * @brief Definition of a class to assemble matrix entries in a distributed way.
 * @author Thomas Brandes
 * @date 07.09.2017
 */

#pragma once 

#include <scai/lama/storage/COOStorage.hpp>
#include <scai/dmemo/Distribution.hpp>
#include <scai/hmemo/HArray.hpp>
#include <scai/common/Printable.hpp>

#include <vector>

namespace scai
{

namespace lama
{

/** This classs allows to assembly matrix entries by different processors. Each processor can 
 *  add matrix entries by global coordinates and the value. 
 *
 *  An object stands for a set of collected matrix entries and can be converted into a matrix
 *  with a given distribution.
 */
template<typename ValueType>
class MatrixAssembly : public common::Printable
{
public:

    /** Construct an empty matrix assembly
     *
     *  @param comm specifies the processor set onto which the data will be assembled
     */
    MatrixAssembly( const dmemo::CommunicatorPtr comm = dmemo::Communicator::getCommunicatorPtr() );

    /** Destructor of the access, also inserts the assembled entries into the matrix. */

    ~MatrixAssembly();

    /** Might be used for optimization to indicate how many elements might be added by this processor. */

    void reserve( const IndexType n );

    /**
     *  @brief Remove all assembled entries.
     *
     *  This method might be helpful to use an object of this class for
     *  multiple assembling phases.
     */
    void clear();

    /** Add a matrix element with global coordinates */

    void push( const IndexType i, const IndexType j, const ValueType val );

    /** Get the number of rows used in the assembled data.
     *
     *  For conversion into a matrix this matrix must have at least the 
     *  size getNumRows() $x$ getNumColumns()
     */

    IndexType getNumRows() const;

    /** Get the number of columns used in the assembled data.
     *
     *  For conversion into a matrix this matrix must have at least the 
     *  size getNumRows() $x$ getNumColumns()
     */

    IndexType getNumColumns() const;

    /** Get the number of entries collected in the assembly so far.  
     *
     */
    IndexType getNumValues() const;

    /** Getter for the communicator used for the assembly matrix. */

    const dmemo::Communicator& getCommunicator() const;

    /**
     *  @brief Get the assembled data in COO format localized for a given distribution
     *
     *  Be careful as the returned coo storage can have multiple entries for the same matrix position.
     */
    COOStorage<ValueType> buildLocalCOO( const dmemo::Distribution& dist, const IndexType numColumns ) const;

    /**
     *  @brief Get the assembled data in COO format replicated on all processors
     *
     *  Be careful as the returned coo storage can have multiple entries for the same matrix position.
     */
    COOStorage<ValueType> buildGlobalCOO( const IndexType numRows, const IndexType numColumns ) const;

    /**
     *  Override Printable::writeAt 
     */
    virtual void writeAt( std::ostream& stream ) const;

private:

    SCAI_LOG_DECL_STATIC_LOGGER( logger )

    /** Resort COO data according to the ownership of the row indexes */

    void exchangeCOO(                      
        hmemo::HArray<IndexType>& outIA,
        hmemo::HArray<IndexType>& outJA,
        hmemo::HArray<ValueType>& outValues,
        const hmemo::HArray<IndexType> inIA,
        const hmemo::HArray<IndexType> inJA,
        const hmemo::HArray<ValueType> inValues,
        const dmemo::Distribution& dist ) const;

    /** Check for correct indexes of the assembled entries */

    void checkLegalIndexes( const IndexType numRows, const IndexType numColumns ) const;

    // for pushing the assembled data we use the C++ vector class

    std::vector<IndexType> mIA;
    std::vector<IndexType> mJA;
    std::vector<ValueType> mValues;

    dmemo::CommunicatorPtr mComm;
};

/* ================================================================================ */
/*   Implementation of inline methods                                               */
/* ================================================================================ */

template<typename ValueType>
void MatrixAssembly<ValueType>::push( const IndexType i, const IndexType j, const ValueType val )
{
    mIA.push_back( i );
    mJA.push_back( j );
    mValues.push_back( val );
}   

template<typename ValueType>
const dmemo::Communicator& MatrixAssembly<ValueType>::getCommunicator() const
{
    return *mComm;
}

}

}
