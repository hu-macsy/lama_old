/**
 * @file CSRStorage.cpp
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
 * @brief Implementation and instantiation for template class CSRStorage.
 * @author Thomas Brandes
 * @date 04.06.2011
 */

// hpp
#include <scai/lama/storage/CSRStorage.hpp>
#include <scai/lama/storage/COOStorage.hpp>

// local library
#include <scai/utilskernel/UtilKernelTrait.hpp>
#include <scai/sparsekernel/CSRKernelTrait.hpp>
#include <scai/sparsekernel/CSRUtils.hpp>
#include <scai/sparsekernel/COOUtils.hpp>
#include <scai/sparsekernel/DIAKernelTrait.hpp>

#include <scai/lama/storage/StorageMethods.hpp>
#include <scai/lama/Scalar.hpp>

#include <scai/dmemo/Redistributor.hpp>


// internal scai libraries
#include <scai/sparsekernel/openmp/OpenMPCSRUtils.hpp>
#include <scai/utilskernel/openmp/OpenMPUtils.hpp>
#include <scai/utilskernel/HArrayUtils.hpp>
#include <scai/utilskernel/LAMAKernel.hpp>

#include <scai/blaskernel/BLASKernelTrait.hpp>

#include <scai/hmemo.hpp>

#include <scai/tracing.hpp>

#include <scai/common/macros/assert.hpp>
#include <scai/common/Constants.hpp>
#include <scai/common/macros/print_string.hpp>
#include <scai/common/macros/unsupported.hpp>
#include <scai/tasking/NoSyncToken.hpp>
#include <scai/common/TypeTraits.hpp>
#include <scai/common/Math.hpp>
#include <scai/common/macros/instantiate.hpp>

#include <memory>

using std::shared_ptr;

namespace scai
{

using namespace hmemo;
using namespace dmemo;
using namespace utilskernel;

using std::unique_ptr;
using common::TypeTraits;
using common::BinaryOp;
using common::CompareOp;

using sparsekernel::CSRKernelTrait;
using sparsekernel::DIAKernelTrait;
using sparsekernel::OpenMPCSRUtils;
using sparsekernel::CSRUtils;

using tasking::SyncToken;

namespace lama
{

/* --------------------------------------------------------------------------- */

SCAI_LOG_DEF_TEMPLATE_LOGGER( template<typename ValueType>, CSRStorage<ValueType>::logger, "MatrixStorage.CSRStorage" )

/* --------------------------------------------------------------------------- */

template<typename ValueType>
CSRStorage<ValueType>::CSRStorage( ContextPtr ctx ) :

    MatrixStorage<ValueType>( 0, 0, ctx ),

    mIA( IndexType( 1 ), IndexType( 0 ), ctx ),
    mJA( ctx ),
    mValues( ctx )
{
    // no row indexes necessary
 
    mDiagonalProperty = checkDiagonalProperty();
    mSortedRows = false;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
CSRStorage<ValueType>::CSRStorage( IndexType numRows, IndexType numColumns, ContextPtr ctx ) :

    MatrixStorage<ValueType>( numRows, numColumns, ctx ),
    mIA( numRows + 1, IndexType( 0 ), ctx ),
    mJA( ctx ),
    mValues( ctx )
{
    SCAI_LOG_DEBUG( logger, "COOStorage for matrix " << getNumRows()
                             << " x " << getNumColumns() << ", no non-zero elements @ " << *ctx )

    mSortedRows = false;

    _MatrixStorage::resetDiagonalProperty();
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
CSRStorage<ValueType>::CSRStorage(
    const IndexType numRows,
    const IndexType numColumns,
    HArray<IndexType> ia,
    HArray<IndexType> ja,
    HArray<ValueType> values,
    ContextPtr ctx ) :

    MatrixStorage<ValueType>( numRows, numColumns, ctx ),
    mIA( std::move( ia ) ),
    mJA( std::move( ja ) ),
    mValues( std::move( values ) )
{
    SCAI_LOG_DEBUG( logger, "input data: ia = " << ia << ", ja = " << ja << ", values = " << values )
    SCAI_LOG_DEBUG( logger, "my data: mIA = " << mIA << ", mJA = " << mJA << ", values = " << mValues )

    SCAI_LOG_INFO( logger, "check valid CSR arrays for constructor" )

    // full consistency checks for valid CSR data

    SCAI_ASSERT_EQ_ERROR( mIA.size(), numRows + 1, "csr offset array has illegal size" )
    SCAI_ASSERT_ERROR( HArrayUtils::isSorted( mIA, CompareOp::LE ), "ia is invalid offset array, entries not ascending" )
    SCAI_ASSERT_EQ_ERROR( mJA.size(), HArrayUtils::getVal( mIA, numRows ), "last entry in offsets must be size of ja" );

    SCAI_ASSERT_EQ_ERROR( mJA.size(), mValues.size(), "serious mismatch of CSR ja and values array" )
    SCAI_ASSERT_ERROR( HArrayUtils::validIndexes( mJA, numColumns ), "invalid column indexes, #cols = " << numColumns );

    // now set properties

    mDiagonalProperty = checkDiagonalProperty();
    mSortedRows       = false;
    buildRowIndexes();
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
CSRStorage<ValueType>& CSRStorage<ValueType>::operator=( const CSRStorage<ValueType>& other )
{
    assignCSR( other );
    return *this;
}

template<typename ValueType>
CSRStorage<ValueType>& CSRStorage<ValueType>::operator=( CSRStorage<ValueType>&& other )
{
    // move of all member variables

    mIA = std::move( other.mIA );
    mJA = std::move( other.mJA );
    mValues = std::move( other.mValues );

    mSortedRows = other.mSortedRows;

    // call of move assignment for base class, use moveImpl

    MatrixStorage<ValueType>::moveImpl( std::move( other ) );

    return *this;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void CSRStorage<ValueType>::splitUp(
    IndexType& numRows,
    IndexType& numColumns,
    hmemo::HArray<IndexType>& ia,
    hmemo::HArray<IndexType>& ja,
    hmemo::HArray<ValueType>& values )
{
    ia = std::move( mIA );
    ja = std::move( mJA );
    values = std::move( mValues );

    // reset the dimensions of this storage to zero so it remains consistent

    _MatrixStorage::splitUp( numRows, numColumns );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void CSRStorage<ValueType>::print( std::ostream& stream ) const
{
    using std::endl;
    stream << "CSRStorage " << getNumRows() << " x " << getNumColumns() << ", #values = " << getNumValues() << endl;

    ContextPtr host = Context::getHostPtr();
    ReadAccess<IndexType> ia( mIA, host );
    ReadAccess<IndexType> ja( mJA, host );
    ReadAccess<ValueType> values( mValues, host );

    for ( IndexType i = 0; i < getNumRows(); i++ )
    {
        stream << "Row " << i << " ( " << ia[i] << " - " << ia[i + 1] << " ) :";

        for ( IndexType jj = ia[i]; jj < ia[i + 1]; ++jj )
        {
            stream << " " << ja[jj] << ":" << values[jj];
        }

        stream << endl;
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void CSRStorage<ValueType>::clear()
{
    allocate( 0, 0 );   // sets everything correctly, ia array has one entry
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
Format CSRStorage<ValueType>::getFormat() const
{
    return Format::CSR;
}

/* --------------------------------------------------------------------------- */

#ifndef SCAI_ASSERT_LEVEL_DEBUG
template<typename ValueType>
void CSRStorage<ValueType>::check( const char* ) const
{}
#else
template<typename ValueType>
void CSRStorage<ValueType>::check( const char* msg ) const
{
    SCAI_ASSERT_EQUAL_ERROR( getNumRows() + 1, mIA.size() )
    SCAI_ASSERT_EQ_ERROR( mJA.size(), mValues.size(), "serious mistmach for sizes of csrJa and csrValues" )
    // check ascending values in offset array mIA
    {
        static LAMAKernel<UtilKernelTrait::isSorted<IndexType> > isSorted;
        static LAMAKernel<UtilKernelTrait::getValue<IndexType> > getValue;
        ContextPtr loc = this->getContextPtr();
        isSorted.getSupportedContext( loc, getValue );
        ReadAccess<IndexType> csrIA( mIA, loc );
        SCAI_CONTEXT_ACCESS( loc )
        IndexType numValues = getValue[ loc ]( csrIA.get(), getNumRows() );
        SCAI_ASSERT_EQ_ERROR( numValues, mJA.size(),
                              "ia[" << getNumRows() << "] = " << numValues << ", msg = " << msg )
        SCAI_ASSERT_ERROR( isSorted[ loc ]( csrIA.get(), getNumRows() + 1, CompareOp::LE ),
                           *this << " @ " << msg << ": IA is illegal offset array" )
    }
    // check column indexes in JA
    {
        static LAMAKernel<UtilKernelTrait::validIndexes> validIndexes;
        ContextPtr loc = this->getContextPtr();
        validIndexes.getSupportedContext( loc );
        ReadAccess<IndexType> rJA( mJA, loc );
        SCAI_CONTEXT_ACCESS( loc )
        SCAI_ASSERT_ERROR( validIndexes[loc]( rJA.get(), mJA.size(), getNumColumns() ),
                           *this << " @ " << msg << ": illegel indexes in JA" )
    }
}
#endif

/* --------------------------------------------------------------------------- */

template<typename ValueType>
bool CSRStorage<ValueType>::checkDiagonalProperty() const
{
    const IndexType m = getNumRows();
    const IndexType n = getNumColumns();

    // diagonal property is given if size of matrix is 0

    if ( m == 0 || n == 0 )
    {
        return true;
    }

    // non-zero sized matrix with no values has not diagonal property

    if ( mJA.size() == 0 )
    {
        return false;
    }

    return CSRUtils::hasDiagonalProperty( m, n, mIA, mJA, mSortedRows, getContextPtr() );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void CSRStorage<ValueType>::setIdentity( const IndexType size )
{
    SCAI_LOG_DEBUG( logger, "set identity, size = " << size )

    _MatrixStorage::setDimension( size, size );

    // Note: we pass also the current context as it might have been changed

    IndexType start = 0;
    IndexType inc   = 1;

    HArrayUtils::setSequence( mIA, start, inc, size + 1, getContextPtr() );
    HArrayUtils::setSequence( mJA, start, inc, size, getContextPtr() );
    HArrayUtils::setSameValue( mValues, size, ValueType( 1 ), getContextPtr() );

    mDiagonalProperty = true;  // obviously given for identity matrix
    mSortedRows = true;        // obviously given for identity matrix

    // Note: we do not build row indexes, no row is empty

    SCAI_LOG_INFO( logger, *this << ": identity matrix" )
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void CSRStorage<ValueType>::assignDiagonal( const HArray<ValueType>& diagonal )
{
    const IndexType size = diagonal.size();

    _MatrixStorage::setDimension( size, size );

    // Note: we pass also the current context as it might have been changed
    
    IndexType start = 0;
    IndexType inc   = 1;

    HArrayUtils::setSequence( mIA, start, inc, size + 1, getContextPtr() );
    HArrayUtils::setSequence( mJA, start, inc, size, getContextPtr() );
    HArrayUtils::setArray( mValues, diagonal, common::BinaryOp::COPY, getContextPtr() );

    mDiagonalProperty = true; // obviously given for identity matrix
    mSortedRows = true; // obviously given for identity matrix

    // Note: we do not build row indexes, no row is empty

    SCAI_LOG_INFO( logger, *this << ": diagonal matrix" )
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
template<typename OtherValueType>
void CSRStorage<ValueType>::setCSRDataImpl(
    const IndexType numRows,
    const IndexType numColumns,
    const HArray<IndexType>& ia,
    const HArray<IndexType>& ja,
    const HArray<OtherValueType>& values,
    const ContextPtr /* loc */ )
{
    SCAI_REGION( "Storage.CSR.setCSR" )

    ContextPtr loc = this->getContextPtr();

    IndexType numValues;

    if ( ia.size() == numRows )
    {
        numValues = HArrayUtils::reduce( ia, BinaryOp::ADD );
    }
    else if ( ia.size() == numRows + 1 )
    {
        SCAI_ASSERT( HArrayUtils::isSorted( ia, CompareOp::LE ), "ia is invalid offset array, entries not ascending" )
        numValues= HArrayUtils::getVal( ia, numRows );  // last entry is num values
    }
    else
    {
        COMMON_THROWEXCEPTION( "ia array with size = " << ia.size() << " illegal, #rows = " << numRows )
    }

    SCAI_ASSERT_EQ_ERROR( numValues, ja.size(), "ja array has illegal size" );
    SCAI_ASSERT_EQ_ERROR( numValues, values.size(), "ja array has illegal size" );

    SCAI_ASSERT( HArrayUtils::validIndexes( ja, numColumns ), "invalid column indexes, #cols = " << numColumns );

    // now we can copy all data

    _MatrixStorage::setDimension( numRows, numColumns );

    SCAI_LOG_DEBUG( logger, "fill " << *this << " with csr data, " << numValues << " non-zero values" )

    // storage data will be directly allocated on the location

    if ( ia.size() == numRows )
    {
        {
            // reserve enough memory for mIA
            WriteOnlyAccess<IndexType> myIA( mIA, loc, numRows + 1 );
        }
        HArrayUtils::assign( mIA, ia, loc );
        HArrayUtils::scan1( mIA, loc );
    }
    else
    {
        HArrayUtils::assign( mIA, ia, loc );
    }

    HArrayUtils::assign( mValues, values, loc );
    HArrayUtils::assign( mJA, ja, loc );

    /* do not sort rows, destroys diagonal property during redistribute

     OpenMPCSRUtils::sortRowElements( myJA.get(), myValues.get(), myIA.get(),
     getNumRows(), mDiagonalProperty );

     */
    mDiagonalProperty = checkDiagonalProperty();
    mSortedRows       = false;
    buildRowIndexes();
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void CSRStorage<ValueType>::sortRows( bool keepDiagonalFirst )
{
    // sort in place 

    CSRUtils::sortRows( mJA, mValues, mIA, getNumRows(), getNumColumns(), keepDiagonalFirst, getContextPtr() );

    if ( keepDiagonalFirst )
    {
        mSortedRows       = false;     // no optimization for routines that rely on sorted data
    }
    else
    {
        mSortedRows       = true;
    }

    // Note: diagonal property remains unchanged
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
bool CSRStorage<ValueType>::hasSortedRows( bool allowDiagonalFirst )
{
    // sort in place 

    return sparsekernel::CSRUtils::hasSortedRows( mIA, mJA, getNumRows(), getNumColumns(), allowDiagonalFirst, getContextPtr() );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void CSRStorage<ValueType>::setDiagonalFirst()
{
    IndexType numDiagonals = std::min( getNumRows(), getNumColumns() );
    IndexType numFirstDiagonals = sparsekernel::CSRUtils::setDiagonalFirst( mJA, mValues, mIA, getNumColumns(), getContextPtr() );

    mDiagonalProperty = numDiagonals == numFirstDiagonals;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void CSRStorage<ValueType>::setDiagonalProperty()
{
    setDiagonalFirst();

    SCAI_ASSERT( mDiagonalProperty, "Missing diagonal element, cannot set diagonal property" )
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void CSRStorage<ValueType>::buildRowIndexes()
{
    mRowIndexes.clear();

    if ( mDiagonalProperty )
    {
        SCAI_LOG_INFO( logger, "buildRowIndexes not done as diagonal property implies at least one entry per row" );
        return;
    }

    if ( getNumRows() == 0 )
    {
        return;
    }

    // Build the row indexes if the ratio of non-zero rows is below the threshold

    sparsekernel::CSRUtils::nonEmptyRows( mRowIndexes, mIA, mCompressThreshold, getContextPtr() );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void CSRStorage<ValueType>::redistributeCSR( const CSRStorage<ValueType>& other, const dmemo::Redistributor& redistributor )
{
    SCAI_REGION( "Storage.redistributeCSR" )

    const dmemo::Distribution& sourceDistribution = *redistributor.getSourceDistributionPtr();
    const dmemo::Distribution& targetDistribution = *redistributor.getTargetDistributionPtr();
    SCAI_LOG_INFO( logger,
                   other << ": redistribute of CSR<" << other.getValueType() << "> to CSR<" << this->getValueType() << " via " << redistributor )
    bool sameDist = false;

    // check for same distribution, either equal or both replicated

    if ( sourceDistribution.isReplicated() && targetDistribution.isReplicated() )
    {
        sameDist = true;
    }
    else if ( &sourceDistribution == &targetDistribution )
    {
        sameDist = true;
    }

    if ( sameDist )
    {
        SCAI_LOG_INFO( logger, "redistributor with same source/target distribution" )
        assign( other );
        return; // so we are done
    }

    // check that source distribution fits with storage
    SCAI_ASSERT_EQUAL_ERROR( other.getNumRows(), sourceDistribution.getLocalSize() )

    if ( &other == this )
    {
        // due to alias we need temporary array
        HArray<IndexType> targetIA;
        HArray<IndexType> targetJA;
        HArray<ValueType> targetValues;
        StorageMethods<ValueType>::redistributeCSR( targetIA, targetJA, targetValues, other.getIA(), other.getJA(),
                other.getValues(), redistributor );
        // we can swap the new arrays
        mIA.swap( targetIA );
        mJA.swap( targetJA );
        mValues.swap( targetValues );
    }
    else
    {
        StorageMethods<ValueType>::redistributeCSR( mIA, mJA, mValues, other.getIA(), other.getJA(), other.getValues(),
                redistributor );
    }

    // it is not necessary to convert the other storage to CSR
 
    _MatrixStorage::setDimension( mIA.size() - 1, other.getNumColumns() );

    mDiagonalProperty = checkDiagonalProperty();
    buildRowIndexes();
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
CSRStorage<ValueType>::~CSRStorage()
{
    SCAI_LOG_DEBUG( logger,
                    "~CSRStorage, size = " << getNumRows() << " x " << getNumColumns() << ", # non-zeros = " << getNumValues() )
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
IndexType CSRStorage<ValueType>::getNumValues() const
{
    return mJA.size();
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void CSRStorage<ValueType>::purge()
{
    // delete all old values
    mIA.purge();
    mJA.purge();
    mValues.purge();
    allocate( 0, 0 );   // sets everything correctly, ia array has one entry
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void CSRStorage<ValueType>::allocate( IndexType numRows, IndexType numColumns )
{
    SCAI_LOG_INFO( logger,
                   "allocate CSR sparse matrix of size " << numRows << " x " << numColumns << ", numValues = 0" )

    _MatrixStorage::setDimension( numRows, numColumns );

    mJA.clear();
    mValues.clear();
    mIA.clear();

    // offset array requires initialization to have consistent data

    mIA.resize( getNumRows() + 1 );
    HArrayUtils::setScalar( mIA, IndexType( 0 ), common::BinaryOp::COPY, getContextPtr() );

    mDiagonalProperty = checkDiagonalProperty();
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void CSRStorage<ValueType>::compress( const RealType<ValueType> eps, bool keepDiagonal )
{
    IndexType numValues = mJA.size();

    sparsekernel::CSRUtils::compress( mIA, mJA, mValues, keepDiagonal, eps, this->getContextPtr() );

    if ( mJA.size() == numValues )
    {
        SCAI_LOG_INFO( logger, "compress, nothing removed" )
    }
    else
    {
        SCAI_LOG_INFO( logger, "compress, now " << mJA.size() << " entries left from " << numValues )

        // sorting has not changed

        mDiagonalProperty = checkDiagonalProperty();
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void CSRStorage<ValueType>::swap( CSRStorage<ValueType>& other )
{
    // swap data of base class

    MatrixStorage<ValueType>::swap( other );

    // swap own member variables

    std::swap( mSortedRows, other.mSortedRows );

    mIA.swap( other.mIA );
    mJA.swap( other.mJA );
    mValues.swap( other.mValues );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void CSRStorage<ValueType>::swap( HArray<IndexType>& ia, HArray<IndexType>& ja, HArray<ValueType>& values )
{
    SCAI_ASSERT_EQ_ERROR( ia.size(), getNumRows() + 1, "swap with illegally sized ia array" )

    IndexType numValues = ia[ getNumRows() ];

    SCAI_ASSERT_EQUAL_ERROR( numValues, ja.size() )
    SCAI_ASSERT_EQUAL_ERROR( numValues, values.size() )
    mIA.swap( ia );
    mJA.swap( ja );
    mValues.swap( values );
    mDiagonalProperty = checkDiagonalProperty();
    // build new array of row indexes
    buildRowIndexes();
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
size_t CSRStorage<ValueType>::getMemoryUsageImpl() const
{
    size_t memoryUsage = 0;
    memoryUsage += sizeof( IndexType );
    memoryUsage += sizeof( IndexType ) * mIA.size();
    memoryUsage += sizeof( IndexType ) * mJA.size();
    memoryUsage += sizeof( ValueType ) * mValues.size();
    return memoryUsage;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void CSRStorage<ValueType>::writeAt( std::ostream& stream ) const
{
    stream << "CSRStorage<" << common::getScalarType<ValueType>() << ">("
           << " size = " << getNumRows() << " x " << getNumColumns()
           << ", nnz = " << getNumValues() << ", diag = " << mDiagonalProperty
           << ", sorted = " << mSortedRows << ", ctx = " << *getContextPtr() << " )";
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
ValueType CSRStorage<ValueType>::getValue( const IndexType i, const IndexType j ) const
{
    SCAI_ASSERT_VALID_INDEX_DEBUG( i, getNumRows(), "row index out of range" )
    SCAI_ASSERT_VALID_INDEX_DEBUG( j, getNumColumns(), "column index out of range" )

    SCAI_LOG_TRACE( logger, "get value (" << i << ", " << j << ")" )

    static LAMAKernel<CSRKernelTrait::getValuePos> getValuePos;

    ContextPtr loc = this->getContextPtr();
    getValuePos.getSupportedContext( loc );
    SCAI_CONTEXT_ACCESS( loc )

    ReadAccess<IndexType> rIa( mIA, loc );
    ReadAccess<IndexType> rJa( mJA, loc );

    IndexType pos = getValuePos[loc]( i, j, rIa.get(), rJa.get() );

    ValueType val = 0;

    if ( pos != invalidIndex )
    {
        SCAI_ASSERT_VALID_INDEX_DEBUG( pos, mJA.size(), "illegal value position for ( " << i << ", " << j << " )" );

        val = mValues[ pos ];
    }

    return val;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void CSRStorage<ValueType>::setValue( const IndexType i,
                                      const IndexType j,
                                      const ValueType val,
                                      const BinaryOp op )
{
    SCAI_ASSERT_VALID_INDEX_DEBUG( i, getNumRows(), "row index out of range" )
    SCAI_ASSERT_VALID_INDEX_DEBUG( j, getNumColumns(), "column index out of range" )

    SCAI_LOG_DEBUG( logger, "set value (" << i << ", " << j << ")" )

    static LAMAKernel<CSRKernelTrait::getValuePos> getValuePos;

    ContextPtr loc = this->getContextPtr();
    getValuePos.getSupportedContext( loc );

    SCAI_CONTEXT_ACCESS( loc )

    ReadAccess<IndexType> rIa( mIA, loc );
    ReadAccess<IndexType> rJa( mJA, loc );

    IndexType pos = getValuePos[loc]( i, j, rIa.get(), rJa.get() );

    if ( pos == invalidIndex )
    {
        COMMON_THROWEXCEPTION( "CSR storage has no entry ( " << i << ", " << j << " ) " )
    }

    utilskernel::HArrayUtils::setVal( mValues, pos, val, op );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void CSRStorage<ValueType>::prefetch( const ContextPtr location ) const
{
    mIA.prefetch( location );
    mJA.prefetch( location );
    mValues.prefetch( location );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
const HArray<IndexType>& CSRStorage<ValueType>::getIA() const
{
    return mIA;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
const HArray<IndexType>& CSRStorage<ValueType>::getJA() const
{
    return mJA;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
const HArray<ValueType>& CSRStorage<ValueType>::getValues() const
{
    return mValues;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void CSRStorage<ValueType>::setDiagonal( const ValueType value )
{
    bool okay = sparsekernel::CSRUtils::setDiagonal( mValues, value, getNumRows(), getNumColumns(), 
                                                     mIA, mJA, mSortedRows, getContextPtr() );

    if ( !okay )
    {
        COMMON_THROWEXCEPTION( "could not set diagonal, not all entries available" )
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void CSRStorage<ValueType>::setDiagonalV( const HArray<ValueType>& diagonal )
{
    bool okay = sparsekernel::CSRUtils::setDiagonalV( mValues, diagonal, getNumRows(), getNumColumns(), 
                                                      mIA, mJA, mSortedRows, getContextPtr() );

    if ( !okay )
    {
        COMMON_THROWEXCEPTION( "could not set diagonal, not all entries available" )
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void CSRStorage<ValueType>::getRow( HArray<ValueType>& row, const IndexType i ) const
{
    SCAI_REGION( "Storage.CSR.getRow" )

    SCAI_ASSERT_VALID_INDEX_DEBUG( i, getNumRows(), "row index out of range" )

    row.setSameValue( getNumColumns(), ValueType( 0 ) );

    IndexType n1    = mIA[i];
    IndexType nrow  = mIA[i + 1] - n1;

    // row [ mJA[n1:] ] = mValues[n1:],

    static LAMAKernel<UtilKernelTrait::setScatter<ValueType, ValueType> > setScatter;

    ContextPtr loc = this->getContextPtr();
    setScatter.getSupportedContext( loc );

    SCAI_CONTEXT_ACCESS( loc )

    WriteAccess<ValueType> wRow( row, loc, getNumColumns() );
    const ReadAccess<IndexType> ja( mJA, loc );
    const ReadAccess<ValueType> values( mValues, loc );

    setScatter[loc]( wRow.get(), ja.get() + n1, true, values.get() + n1, BinaryOp::COPY, nrow );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void CSRStorage<ValueType>::getSparseRow( hmemo::HArray<IndexType>& jA, hmemo::HArray<ValueType>& values, const IndexType i ) const
{
    SCAI_REGION( "Storage.CSR.getSparseRow" )

    SCAI_ASSERT_VALID_INDEX_DEBUG( i, getNumRows(), "row index out of range" )

    IndexType offs  = mIA[i];
    IndexType nrow  = mIA[i + 1] - offs;

    // resize the output arrays, invalidate old data before

    jA.clear();
    jA.resize( nrow );
    values.clear();
    values.resize( nrow );

    // just copy the corresponding parts of the csrJA and csrValues array

    BinaryOp op = BinaryOp::COPY;

    const IndexType inc = 1;  // no strides

    HArrayUtils::setArraySection( jA, 0, inc, mJA, offs, inc, nrow, op, getContextPtr() );
    HArrayUtils::setArraySection( values, 0, inc, mValues, offs, inc, nrow, op, getContextPtr() );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void CSRStorage<ValueType>::getSparseColumn( 
    hmemo::HArray<IndexType>& iA, 
    hmemo::HArray<ValueType>& values, 
    const IndexType j ) const
{
    SCAI_REGION( "Storage.CSR.getSparseCol" )

    SCAI_ASSERT_VALID_INDEX_DEBUG( j, getNumColumns(), "col index out of range" )

    static LAMAKernel<CSRKernelTrait::getValuePosCol> getValuePosCol;

    ContextPtr loc = this->getContextPtr();

    getValuePosCol.getSupportedContext( loc );

    HArray<IndexType> valuePos;     // positions in the values array

    {
        SCAI_CONTEXT_ACCESS( loc )

        WriteOnlyAccess<IndexType> wRowIndexes( iA, loc, getNumRows() );
        WriteOnlyAccess<IndexType> wValuePos( valuePos, loc, getNumRows() );

        ReadAccess<IndexType> rIA( mIA, loc );
        ReadAccess<IndexType> rJA( mJA, loc );

        IndexType cnt = getValuePosCol[loc]( wRowIndexes.get(), wValuePos.get(), j,
                                             rIA.get(), getNumRows(), rJA.get(), rJA.size() );

        wRowIndexes.resize( cnt );
        wValuePos.resize( cnt );
    }

    // values = mValues[ pos ];

    HArrayUtils::gather( values, mValues, valuePos, BinaryOp::COPY, loc );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void CSRStorage<ValueType>::setRow( const HArray<ValueType>& row, const IndexType i, const BinaryOp op )
{
    SCAI_REGION( "Storage.CSR.setRow" )

    SCAI_ASSERT_VALID_INDEX_DEBUG( i, getNumRows(), "row index out of range" )
    SCAI_ASSERT_GE_DEBUG( row.size(), getNumColumns(), "row array to small for set" )

    IndexType n1    = mIA[i];
    IndexType nrow  = mIA[i + 1] - n1;

    // mValues[n1:] op= row [ mJA[n1:] ]

    HArray<ValueType> gatheredRow;  // temporary required as gather does not support op

    ContextPtr loc = this->getContextPtr();

    {
        static LAMAKernel<UtilKernelTrait::setGather<ValueType, ValueType> > setGather;

        setGather.getSupportedContext( loc );

        SCAI_CONTEXT_ACCESS( loc )

        WriteOnlyAccess<ValueType> wRow( gatheredRow, loc, nrow );
        const ReadAccess<IndexType> ja( mJA, loc );
        const ReadAccess<ValueType> rRow( row, loc );

        setGather[loc]( wRow.get(), rRow.get(), ja.get() + n1, BinaryOp::COPY, nrow );
    }

    HArrayUtils::setArraySection( mValues, n1, 1, gatheredRow, 0, 1, nrow, op, loc );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void CSRStorage<ValueType>::getColumn( HArray<ValueType>& column, const IndexType j ) const
{
    SCAI_REGION( "Storage.CSR.getDenseCol" )

    HArray<IndexType> rowIndexes;   // row indexes that have entry for column j
    HArray<ValueType> colValues;    // contains the values of entries belonging to column j

    getSparseColumn( rowIndexes, colValues, j );

    HArrayUtils::buildDenseArray( column, getNumRows(), colValues, rowIndexes, ValueType( 0 ) );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void CSRStorage<ValueType>::setColumn( 
        const HArray<ValueType>& column, 
        const IndexType j,
        const BinaryOp op )
{
    SCAI_REGION( "Storage.CSR.setCol" )
    SCAI_ASSERT_VALID_INDEX_DEBUG( j, getNumColumns(), "column index out of range" )
    SCAI_ASSERT_GE_DEBUG( column.size(), getNumRows(), "column array to small for set" )

    static LAMAKernel<CSRKernelTrait::getValuePosCol> getValuePosCol;

    ContextPtr loc = this->getContextPtr();

    getValuePosCol.getSupportedContext( loc );

    HArray<IndexType> rowIndexes;   // row indexes that have entry for column j
    HArray<IndexType> valuePos;     // positions in the values array
    HArray<ValueType> colValues;    // contains the values of entries belonging to column j

    {
        SCAI_CONTEXT_ACCESS( loc )

        // allocate rowIndexes, valuePos with maximal possible size

        WriteOnlyAccess<IndexType> wRowIndexes( rowIndexes, loc, getNumRows() );
        WriteOnlyAccess<IndexType> wValuePos( valuePos, loc, getNumRows() );
        ReadAccess<IndexType> rIA( mIA, loc );
        ReadAccess<IndexType> rJA( mJA, loc );

        IndexType cnt = getValuePosCol[loc]( wRowIndexes.get(), wValuePos.get(), j,
                                             rIA.get(), getNumRows(), rJA.get(), rJA.size() );

        wRowIndexes.resize( cnt );
        wValuePos.resize( cnt );
    }

    //  mValues[ pos ] op= column[row]

    HArrayUtils::gather( colValues, column, rowIndexes, BinaryOp::COPY, loc );
    HArrayUtils::scatter( mValues, valuePos, true, colValues, op, loc );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void CSRStorage<ValueType>::getDiagonal( HArray<ValueType>& diagonal ) const
{
    sparsekernel::CSRUtils::getDiagonal( diagonal, getNumRows(), getNumColumns(),
                                         mIA, mJA, mValues, mSortedRows, getContextPtr() );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void CSRStorage<ValueType>::scale( const ValueType value )
{
    HArrayUtils::compute( mValues, mValues, BinaryOp::MULT, value, this->getContextPtr() );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void CSRStorage<ValueType>::conj()
{
    HArrayUtils::unaryOp( mValues, mValues, common::UnaryOp::CONJ, this->getContextPtr() );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void CSRStorage<ValueType>::scaleRows( const HArray<ValueType>& diagonal )
{
    IndexType n = std::min( getNumRows(), diagonal.size() );
    static LAMAKernel<CSRKernelTrait::scaleRows<ValueType> > scaleRows;
    ContextPtr loc = this->getContextPtr();
    scaleRows.getSupportedContext( loc );
    SCAI_CONTEXT_ACCESS( loc )
    {
        ReadAccess<ValueType> rDiagonal( diagonal, loc );
        ReadAccess<IndexType> csrIA( mIA, loc );
        WriteAccess<ValueType> csrValues( mValues, loc ); // updateAccess
        scaleRows[loc]( csrValues.get(), csrIA.get(), n, rDiagonal.get() );
    }

    if ( SCAI_LOG_TRACE_ON( logger ) )
    {
        SCAI_LOG_TRACE( logger, "CSR after scale diagonal" )
        print( std::cout );
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void CSRStorage<ValueType>::wait() const
{
    mIA.wait();
    mJA.wait();
    mValues.wait();
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void CSRStorage<ValueType>::assign( const _MatrixStorage& other )
{
    // translate virtual call to specific template call via wrapper

    mepr::StorageWrapper<CSRStorage, SCAI_NUMERIC_TYPES_HOST_LIST>::assignImpl( this, other );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
template<typename OtherValueType>
void CSRStorage<ValueType>::assignImpl( const MatrixStorage<OtherValueType>& other )
{
    if ( static_cast<const _MatrixStorage*>( &other ) == this )
    {
        SCAI_LOG_INFO( logger, typeName() << ": self assign, skipped, matrix = " << other )
    }
    else if ( other.getFormat() == Format::CSR )
    {
        assignCSR( static_cast<const CSRStorage<OtherValueType>&>( other ) );
    }
    else 
    {
        ContextPtr ctx = getContextPtr();   // will force a valid copy in this context

        other.buildCSRData( mIA, mJA, mValues );

        // setCSRDataImpl can deal with alias and takes advantage of it 
        // and it will set all relevant attributes of this storage correctly.

        setCSRDataImpl( other.getNumRows(), other.getNumColumns(), mIA, mJA, mValues, ctx );
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
template<typename OtherValueType>
void CSRStorage<ValueType>::assignCSR( const CSRStorage<OtherValueType>& other )
{
    ContextPtr ctx = getContextPtr();

    // Note: attributes mContextPtr, mEps remain unchanged

    _MatrixStorage::setDimension( other.getNumRows(), other.getNumColumns() );

    HArrayUtils::assign( mIA, other.getIA(), ctx );
    HArrayUtils::assign( mJA, other.getJA(), ctx );
    HArrayUtils::assign( mValues, other.getValues(), ctx );

    mSortedRows = false; // other.mSortedRows;

    _MatrixStorage::resetDiagonalProperty();

    buildRowIndexes();
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void CSRStorage<ValueType>::assignTranspose( const MatrixStorage<ValueType>& other )
{
    SCAI_REGION( "Storage.CSR.assignTranspose" )

    SCAI_LOG_INFO( logger, *this << ": (CSR) assign transpose " << other )

    // pass HArrays of this storage to build the values in it

    if ( &other == this )
    {
        SCAI_LOG_INFO( logger, *this << ": (CSR) assign transpose in place" )
        HArray<IndexType> tmpIA;
        HArray<IndexType> tmpJA;
        HArray<ValueType> tmpValues;
        // do not change sizes before building CSC data
        other.buildCSCData( tmpIA, tmpJA, tmpValues );
        // sizes must be set correctly BEFORE swap
        _MatrixStorage::_assignTranspose( other );
        swap( tmpIA, tmpJA, tmpValues );  // sets all other data correctly
    }
    else
    {
        _MatrixStorage::_assignTranspose( other );
        SCAI_LOG_INFO( logger, *this << ": (CSR) assign transpose " << other )
        other.buildCSCData( mIA, mJA, mValues );
        mDiagonalProperty = checkDiagonalProperty();
        buildRowIndexes();
    }

    // actualize my member variables (class CSRStorage)
    check( "assignTranspose" );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void CSRStorage<ValueType>::copyTo( _MatrixStorage& other ) const
{
    // Compressed sparse column data can be set directly to other matrix
    other.setCSRData( getNumRows(), getNumColumns(), mIA, mJA, mValues );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void CSRStorage<ValueType>::copyBlockTo( _MatrixStorage& other, const IndexType first, const IndexType n ) const
{
    using namespace utilskernel;

    SCAI_ASSERT_LE( first, first + n, "illegal range" )
    SCAI_ASSERT_VALID_INDEX( first, getNumRows(), "first row out of range" )
    SCAI_ASSERT_VALID_INDEX( first + n - 1, getNumRows(), "last row out of range" );

    // we have not to extract csrIA, csrJA, csrValues, as already available

    ContextPtr loc = this->getContextPtr();

    SCAI_LOG_INFO( logger, "copyBlockTo : first = " << first << ", n = " << n << ", from this : " << *this )

    // copy out the corresponding sections, ia needs a shifting to zero

    HArray<IndexType> blockIA( n + 1 );
    HArrayUtils::setArraySection( blockIA, 0, 1, mIA, first, 1, n +  1, BinaryOp::COPY, loc );

    // Note: BinaryOp::ADD with -offset is strange if IndexType is unsigned

    IndexType offset = blockIA[0];  // gives shifting, as blockIA[0] must be 0
    HArrayUtils::compute( blockIA, blockIA, BinaryOp::SUB, offset, loc );

    IndexType numBlockValues = blockIA[n];

    SCAI_LOG_DEBUG( logger, "offset = " << offset << ", #nnz = " << numBlockValues );

    HArray<IndexType> blockJA( numBlockValues );
    HArray<ValueType> blockValues( numBlockValues );

    HArrayUtils::setArraySection( blockJA, 0, 1, mJA, offset, 1, numBlockValues, BinaryOp::COPY, loc );
    HArrayUtils::setArraySection( blockValues, 0, 1, mValues, offset, 1, numBlockValues, BinaryOp::COPY, loc );

    other.setCSRData( n, getNumColumns(), blockIA, blockJA, blockValues );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
CSRStorage<ValueType>* CSRStorage<ValueType>::newMatrixStorage( const IndexType numRows, const IndexType numColumns ) const
{
    std::unique_ptr<CSRStorage<ValueType> > storage( new CSRStorage<ValueType>( getContextPtr() ) );
    storage->allocate( numRows, numColumns );
    return storage.release();
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void CSRStorage<ValueType>::buildCSRSizes( hmemo::HArray<IndexType>& ia ) const
{
    // build size array from offset array mIA

    static LAMAKernel<CSRKernelTrait::offsets2sizes> offsets2sizes;
    ContextPtr loc = mIA.getValidContext();
    offsets2sizes.getSupportedContext( loc );
    ReadAccess<IndexType> offsetIA( mIA, loc );
    WriteOnlyAccess<IndexType> sizesIA( ia, loc, getNumRows() );
    SCAI_CONTEXT_ACCESS( loc )
    offsets2sizes[ loc ]( sizesIA.get(), offsetIA.get(), getNumRows() );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void CSRStorage<ValueType>::buildCSRData( 
    hmemo::HArray<IndexType>& csrIA, 
    hmemo::HArray<IndexType>& csrJA, 
    hmemo::_HArray& csrValues ) const
{
    mepr::StorageWrapper<CSRStorage, SCAI_NUMERIC_TYPES_HOST_LIST>::
        buildCSRDataImpl( this, csrIA, csrJA, csrValues, getContextPtr() );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
template<typename OtherValueType>
void CSRStorage<ValueType>::buildCSR(
    HArray<IndexType>& ia,
    HArray<IndexType>* ja,
    HArray<OtherValueType>* values,
    const ContextPtr prefLoc ) const
{
    SCAI_REGION( "Storage.CSR.setCSR" )

    // copy the offset array ia and ja

    HArrayUtils::assign( ia, mIA, prefLoc );

    if ( ja != NULL )
    {
        HArrayUtils::assign( *ja, mJA, prefLoc );
    }

    if ( values != NULL )
    {
        HArrayUtils::assign( *values, mValues, prefLoc );
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void CSRStorage<ValueType>::buildCSCData(
    HArray<IndexType>& colIA,
    HArray<IndexType>& colJA,
    HArray<ValueType>& colValues ) const
{
    SCAI_LOG_INFO( logger, *this << ": buildCSCData by call of CSR2CSC" )
    // build the CSC data directly on the device where this matrix is located.
    sparsekernel::CSRUtils::convertCSR2CSC( colIA, colJA, colValues, getNumColumns(), mIA, mJA, mValues, getContextPtr() );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void CSRStorage<ValueType>::splitHalo(
    MatrixStorage<ValueType>& localData,
    MatrixStorage<ValueType>& haloData,
    dmemo::Halo& halo,
    const dmemo::Distribution& colDist,
    const dmemo::Distribution* rowDist ) const
{
    SCAI_REGION( "Storage.splitHalo" )
    SCAI_LOG_INFO( logger, *this << ": split CSR according to column distribution " << colDist )
    SCAI_ASSERT_EQUAL_ERROR( getNumColumns(), colDist.getGlobalSize() )

    if ( colDist.isReplicated() )
    {
        // if there is no column distribution, halo is not needed
        if ( rowDist )
        {
            localData.localize( *this, *rowDist );
        }
        else
        {
            localData.assign( *this );
        }

        haloData.allocate( getNumRows(), 0 );
        halo = Halo(); // empty halo schedule
        return;
    }

    IndexType numRows = getNumRows();

    // check optional row distribution if specified

    if ( rowDist )
    {
        SCAI_LOG_INFO( logger, *this << ": split also localizes for " << *rowDist )
        SCAI_ASSERT_EQUAL_ERROR( getNumRows(), rowDist->getGlobalSize() )
        numRows = rowDist->getLocalSize();
    }

    HArray<IndexType> localIA;
    HArray<IndexType> localJA;
    HArray<ValueType> localValues;
    HArray<IndexType> haloIA;
    HArray<IndexType> haloJA;
    HArray<ValueType> haloValues;
    StorageMethods<ValueType>::splitCSR( localIA, localJA, localValues, haloIA, haloJA, haloValues, mIA, mJA, mValues,
                                         colDist, rowDist );
    SCAI_ASSERT_EQ_DEBUG( localIA.size(), numRows + 1, "illegal size of IA array for local part" )
    SCAI_ASSERT_EQ_DEBUG( haloIA.size(), numRows + 1, "illegal size of IA array for halo part" )
    const IndexType haloNumValues = haloJA.size();
    const IndexType localNumValues = localJA.size();
    SCAI_LOG_INFO( logger,
                   *this << ": split into " << localNumValues << " local non-zeros " " and " << haloNumValues << " halo non-zeros" )
    const IndexType localNumColumns = colDist.getLocalSize();
    IndexType haloNumColumns; // will be available after remap
    // build the halo by the non-local indexes
    _StorageMethods::buildHalo( halo, haloJA, haloNumColumns, colDist );
    SCAI_LOG_INFO( logger, "build halo: " << halo )
    localData.setCSRData( numRows, localNumColumns, localIA, localJA, localValues );
    localData.check( "local part after split" );
    // halo data is expected to have many empty rows, so enable compressing with row indexes
    haloData.setCompressThreshold( 0.5 );
    haloData.setCSRData( numRows, haloNumColumns, haloIA, haloJA, haloValues );
    haloData.check( "halo part after split" );
    SCAI_LOG_INFO( logger,
                   "Result of split: local storage = " << localData << ", halo storage = " << haloData << ", halo = " << halo )
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void CSRStorage<ValueType>::globalizeHaloIndexes( const dmemo::Halo& halo, const IndexType globalNumColumns )
{
    halo.halo2Global( mJA );
    _MatrixStorage::setDimension( getNumRows(), globalNumColumns );
    _MatrixStorage::resetDiagonalProperty();
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void CSRStorage<ValueType>::matrixTimesVector(
    HArray<ValueType>& result,
    const ValueType alpha,
    const HArray<ValueType>& x,
    const ValueType beta,
    const HArray<ValueType>& y,
    const common::MatrixOp op ) const
{
    bool async = false; // synchronously execution, no SyncToken required
    SyncToken* token = gemv( result, alpha, x, beta, y, op, async );
    SCAI_ASSERT( token == NULL, "There should be no sync token for synchronous execution" )
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void CSRStorage<ValueType>::matrixTimesVectorN(
    HArray<ValueType>& result,
    const IndexType n,
    const ValueType alpha,
    const HArray<ValueType>& x,
    const ValueType beta,
    const HArray<ValueType>& y ) const
{
    SCAI_LOG_INFO( logger,
                   *this << ": matrixTimesVector, result = " << result << ", alpha = " << alpha << ", x = " << x << ", beta = " << beta << ", y = " << y )
    SCAI_REGION( "Storage.CSR.timesVectorN" )
    SCAI_ASSERT_EQUAL_ERROR( x.size(), n * getNumColumns() )
    SCAI_ASSERT_EQUAL_ERROR( result.size(), n * getNumRows() )

    if ( ( beta != common::Constants::ZERO ) && ( &result != &y ) )
    {
        SCAI_ASSERT_EQUAL_ERROR( y.size(), n * getNumRows() )
    }

    static LAMAKernel<CSRKernelTrait::gemm<ValueType> > gemm;
    ContextPtr loc = this->getContextPtr();
    gemm.getSupportedContext( loc );
    SCAI_LOG_INFO( logger, *this << ": matrixTimesVectorN on " << *loc )
    ReadAccess<IndexType> csrIA( mIA, loc );
    ReadAccess<IndexType> csrJA( mJA, loc );
    ReadAccess<ValueType> csrValues( mValues, loc );
    ReadAccess<ValueType> rX( x, loc );
    ReadAccess<ValueType> rY( y, loc );
    // due to possible alias of result and y, write access must follow read(y)
    WriteOnlyAccess<ValueType> wResult( result, loc, getNumRows() );
    SCAI_CONTEXT_ACCESS( loc )
    gemm[loc]( wResult.get(), alpha, rX.get(), beta, rY.get(), getNumRows(), n, getNumColumns(),
               csrIA.get(), csrJA.get(), csrValues.get() );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
SyncToken* CSRStorage<ValueType>::sparseGEMV(
    HArray<ValueType>& result,
    const ValueType alpha,
    const HArray<ValueType>& x,
    const common::MatrixOp op,
    bool async ) const
{
    static LAMAKernel<CSRKernelTrait::sparseGEMV<ValueType> > sparseGEMV;
    ContextPtr loc = this->getContextPtr();
    sparseGEMV.getSupportedContext( loc );
    unique_ptr<SyncToken> syncToken;

    if ( async )
    {
        syncToken.reset( loc->getSyncToken() );
    }

    SCAI_ASYNCHRONOUS( syncToken.get() );
    SCAI_CONTEXT_ACCESS( loc )
    ReadAccess<IndexType> csrIA( mIA, loc );
    ReadAccess<IndexType> csrJA( mJA, loc );
    ReadAccess<ValueType> csrValues( mValues, loc );
    ReadAccess<ValueType> rX( x, loc );
    WriteAccess<ValueType> wResult( result, loc );
    // result += alpha * thisMatrix * x, can take advantage of row indexes
    IndexType numNonZeroRows = mRowIndexes.size();
    ReadAccess<IndexType> rows( mRowIndexes, loc );
    sparseGEMV[loc]( wResult.get(), alpha, rX.get(), 
                     numNonZeroRows, rows.get(), 
                     csrIA.get(), csrJA.get(), csrValues.get(), op );

    if ( async )
    {
        syncToken->pushRoutine( rows.releaseDelayed() );
        syncToken->pushRoutine( wResult.releaseDelayed() );
        syncToken->pushRoutine( rX.releaseDelayed() );
        syncToken->pushRoutine( csrIA.releaseDelayed() );
        syncToken->pushRoutine( csrJA.releaseDelayed() );
        syncToken->pushRoutine( csrValues.releaseDelayed() );
    }

    return syncToken.release();
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
SyncToken* CSRStorage<ValueType>::normalGEMV(
    HArray<ValueType>& result,
    const ValueType alpha,
    const HArray<ValueType>& x,
    const ValueType beta,
    const HArray<ValueType>& y,
    const common::MatrixOp op,
    bool async ) const
{
    SCAI_LOG_INFO( logger, "normalGEMV<" << getValueType() << ">( op = " << op << ", async = " << async
                   << ") : result = " << alpha << " * this * x + " << beta << " * y" )

    static LAMAKernel<CSRKernelTrait::normalGEMV<ValueType> > normalGEMV;

    ContextPtr loc = this->getContextPtr();
    normalGEMV.getSupportedContext( loc );
    unique_ptr<SyncToken> syncToken;

    if ( async )
    {
        syncToken.reset( loc->getSyncToken() );
    }

    SCAI_ASYNCHRONOUS( syncToken.get() );
    SCAI_CONTEXT_ACCESS( loc )
    // Note: alias &result == &y possible
    //       ReadAccess on y before WriteOnlyAccess on result guarantees valid data
    ReadAccess<IndexType> csrIA( mIA, loc );
    ReadAccess<IndexType> csrJA( mJA, loc );
    ReadAccess<ValueType> csrValues( mValues, loc );
    ReadAccess<ValueType> rX( x, loc );
    ReadAccess<ValueType> rY( y, loc );
    const IndexType n = common::isTranspose( op ) ? getNumColumns() : getNumRows();
    WriteOnlyAccess<ValueType> wResult( result, loc, n );
    normalGEMV[loc]( wResult.get(), alpha, rX.get(), beta, rY.get(), getNumRows(), getNumColumns(), mJA.size(), csrIA.get(),
                     csrJA.get(), csrValues.get(), op );

    if ( async )
    {
        syncToken->pushRoutine( wResult.releaseDelayed() );
        syncToken->pushRoutine( rY.releaseDelayed() );
        syncToken->pushRoutine( rX.releaseDelayed() );
        syncToken->pushRoutine( csrIA.releaseDelayed() );
        syncToken->pushRoutine( csrJA.releaseDelayed() );
        syncToken->pushRoutine( csrValues.releaseDelayed() );
    }

    return syncToken.release();
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
SyncToken* CSRStorage<ValueType>::normalGEMV(
    HArray<ValueType>& result,
    const ValueType alpha,
    const HArray<ValueType>& x,
    const common::MatrixOp op,
    bool async ) const
{
    SCAI_LOG_INFO( logger, "normalGEMV<" << getValueType() << ">( op = " << op << ", async = " << async
                   << ") : result = " << alpha << " * this * x" )

    static LAMAKernel<CSRKernelTrait::normalGEMV<ValueType> > normalGEMV;
    ContextPtr loc = this->getContextPtr();
    normalGEMV.getSupportedContext( loc );
    unique_ptr<SyncToken> syncToken;

    if ( async )
    {
        syncToken.reset( loc->getSyncToken() );
    }

    SCAI_ASYNCHRONOUS( syncToken.get() );
    SCAI_CONTEXT_ACCESS( loc )
    ReadAccess<IndexType> csrIA( mIA, loc );
    ReadAccess<IndexType> csrJA( mJA, loc );
    ReadAccess<ValueType> csrValues( mValues, loc );
    ReadAccess<ValueType> rX( x, loc );

    const IndexType nTarget = common::isTranspose( op ) ? getNumColumns() : getNumRows();
    WriteOnlyAccess<ValueType> wResult( result, loc, nTarget );
    ValueType beta = 0;
    normalGEMV[loc]( wResult.get(), alpha, rX.get(), beta, NULL, getNumRows(), getNumColumns(), mJA.size(),
                     csrIA.get(), csrJA.get(), csrValues.get(), op );

    if ( async )
    {
        syncToken->pushRoutine( wResult.releaseDelayed() );
        syncToken->pushRoutine( rX.releaseDelayed() );
        syncToken->pushRoutine( csrIA.releaseDelayed() );
        syncToken->pushRoutine( csrJA.releaseDelayed() );
        syncToken->pushRoutine( csrValues.releaseDelayed() );
    }

    return syncToken.release();
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
SyncToken* CSRStorage<ValueType>::gemv(
    HArray<ValueType>& result,
    const ValueType alpha,
    const HArray<ValueType>& x,
    const ValueType beta,
    const HArray<ValueType>& y,
    const common::MatrixOp op,
    bool async ) const
{
    SCAI_REGION( "Storage.CSR.gemv" )

    const IndexType nSource = common::isTranspose( op ) ? getNumRows() : getNumColumns();
    const IndexType nTarget = common::isTranspose( op ) ? getNumColumns() : getNumRows();

    SCAI_LOG_INFO( logger,
                   "GEMV ( op = " << op << ", async = " << async 
                   << " ), result = " << alpha << " * A * x + " << beta << " * y "
                   << ", result = " << result << ", x = " << x << ", y = " << y
                   << ", A (this) = " << *this );

    if ( alpha == common::Constants::ZERO || ( mJA.size() == 0 ) )
    {
        // so we just have result = beta * y, will be done synchronously

        if ( beta == common::Constants::ZERO )
        {
            HArrayUtils::setSameValue( result, nTarget, ValueType( 0 ), this->getContextPtr() );
        }
        else
        {
            HArrayUtils::compute( result, beta, BinaryOp::MULT, y, this->getContextPtr() );
        }

        if ( async )
        {
            return new tasking::NoSyncToken();
        }
        else
        {
            return NULL;
        }
    }

    // check for correct sizes of x

    SCAI_ASSERT_EQ_ERROR( x.size(), nSource, "serious size mismatch gemv, op = " << op )

    if ( beta == common::Constants::ZERO )
    {
        // take version that does not access y at all (can be undefined or aliased to result)

        return normalGEMV( result, alpha, x, op, async );
    }

    // y is relevant, so it must have correct size

    SCAI_ASSERT_EQ_ERROR( y.size(), nTarget, "serious size mismatch gemv, op = " << op )

    if ( &result == &y && ( beta == common::Constants::ONE ) && ( mRowIndexes.size() > 0 ) )
    {
        // y += A * x,  where only some rows in A are filled, uses more efficient routine

        return sparseGEMV( result, alpha, x, op, async );
    }
    else
    {
        return normalGEMV( result, alpha, x, beta, y, op, async );
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
SyncToken* CSRStorage<ValueType>::matrixTimesVectorAsync(
    HArray<ValueType>& result,
    const ValueType alpha,
    const HArray<ValueType>& x,
    const ValueType beta,
    const HArray<ValueType>& y,
    const common::MatrixOp op ) const
{
    bool async = true;
    SyncToken* token = gemv( result, alpha, x, beta, y, op, async );
    SCAI_ASSERT( token, "NULL token not allowed for asynchronous execution gemv, alpha = " << alpha << ", beta = " << beta )
    return token;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void CSRStorage<ValueType>::jacobiIterate(
    HArray<ValueType>& solution,
    const HArray<ValueType>& oldSolution,
    const HArray<ValueType>& rhs,
    const ValueType omega ) const
{
    SCAI_REGION( "Storage.CSR.jacobiIterate" )
    SCAI_LOG_INFO( logger, *this << ": Jacobi iteration for local matrix data." )
    SCAI_ASSERT_ERROR( mDiagonalProperty, *this << ": jacobiIterate requires diagonal property" )

    if ( &solution == &oldSolution )
    {
        COMMON_THROWEXCEPTION( "alias of solution and oldSolution unsupported" )
    }

    SCAI_ASSERT_EQ_DEBUG( getNumRows(), oldSolution.size(), "array oldSolution has illegal size" )
    SCAI_ASSERT_EQ_DEBUG( getNumRows(), rhs.size(), "array rhs has illegal size" )
    SCAI_ASSERT_EQ_DEBUG( getNumRows(), getNumColumns(), "this storage must be square" )

    // matrix must be square
    static LAMAKernel<CSRKernelTrait::jacobi<ValueType> > jacobi;
    ContextPtr loc = this->getContextPtr();
    jacobi.getSupportedContext( loc );
    ReadAccess<IndexType> csrIA( mIA, loc );
    ReadAccess<IndexType> csrJA( mJA, loc );
    ReadAccess<ValueType> csrValues( mValues, loc );
    ReadAccess<ValueType> rOldSolution( oldSolution, loc );
    ReadAccess<ValueType> rRhs( rhs, loc );
    WriteOnlyAccess<ValueType> wSolution( solution, loc, getNumRows() );
    // Due to diagonal property there is no advantage by taking row indexes
    SCAI_CONTEXT_ACCESS( loc )
    jacobi[loc]( wSolution.get(), csrIA.get(), csrJA.get(), csrValues.get(),
                 rOldSolution.get(), rRhs.get(), omega, getNumRows() );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void CSRStorage<ValueType>::jacobiIterateHalo(
    HArray<ValueType>& localSolution,
    const HArray<ValueType>& localDiagonal,
    const HArray<ValueType>& oldHaloSolution,
    const ValueType omega ) const
{
    SCAI_REGION( "Storage.CSR.jacobiIterateHalo" )

    SCAI_LOG_INFO( logger, *this << ": Jacobi iteration for halo matrix data." )

    SCAI_ASSERT_EQ_ERROR( getNumRows(), localSolution.size(), "array localSolution has illegal size" )
    SCAI_ASSERT_EQ_ERROR( getNumColumns(), oldHaloSolution.size(), "array old halo solution has illegal size" )

    static LAMAKernel<CSRKernelTrait::jacobiHalo<ValueType> > jacobiHalo;
    ContextPtr loc = this->getContextPtr();
    jacobiHalo.getSupportedContext( loc );
    {
        WriteAccess<ValueType> wSolution( localSolution, loc ); // will be updated
        ReadAccess<ValueType> localDiagValues( localDiagonal, loc );
        ReadAccess<IndexType> haloIA( mIA, loc );
        ReadAccess<IndexType> haloJA( mJA, loc );
        ReadAccess<ValueType> haloValues( mValues, loc );
        ReadAccess<ValueType> rOldHaloSolution( oldHaloSolution, loc );
        const IndexType numNonEmptyRows = mRowIndexes.size();
        SCAI_LOG_INFO( logger, "#row indexes = " << numNonEmptyRows )

        if ( numNonEmptyRows != 0 )
        {
            ReadAccess<IndexType> haloRowIndexes( mRowIndexes, loc );
            SCAI_CONTEXT_ACCESS( loc )
            jacobiHalo[loc]( wSolution.get(), localDiagValues.get(), haloIA.get(), haloJA.get(), haloValues.get(),
                             haloRowIndexes.get(), rOldHaloSolution.get(), omega, numNonEmptyRows );
        }
        else
        {
            SCAI_CONTEXT_ACCESS( loc )
            jacobiHalo[loc]( wSolution.get(), localDiagValues.get(), haloIA.get(), haloJA.get(), haloValues.get(),
                             NULL, rOldHaloSolution.get(), omega, getNumRows() );
        }
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void CSRStorage<ValueType>::matrixPlusMatrix(
    const ValueType alpha,
    const MatrixStorage<ValueType>& a,
    const ValueType beta,
    const MatrixStorage<ValueType>& b )
{
    SCAI_LOG_INFO( logger, "this = " << alpha << " * A + " << beta << " * B" << ", with A = " << a << ", B = " << b )

    SCAI_REGION( "Storage.CSR.plusMatrix" )

    if ( a.getFormat() != Format::CSR )
    {   
        matrixPlusMatrix( alpha, convert<CSRStorage<ValueType>>( a, getContextPtr() ), beta, b );
    } 
    else if ( b.getFormat() != Format::CSR )
    {   
        matrixPlusMatrix( alpha, a, beta, convert<CSRStorage<ValueType>>( b, getContextPtr() ) );
    }
    else
    {   
        // a and b have the right format so we just cast it correctly

        matrixPlusMatrixImpl( alpha, static_cast<const CSRStorage<ValueType>&>( a ),
                              beta, static_cast<const CSRStorage<ValueType>&>( b ) );
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void CSRStorage<ValueType>::matrixTimesMatrix(
    const ValueType alpha,
    const MatrixStorage<ValueType>& a,
    const MatrixStorage<ValueType>& b,
    const ValueType beta,
    const MatrixStorage<ValueType>& c )
{
    SCAI_LOG_INFO( logger,
                   "this = " << alpha << " +* A * B + " << beta << " * C, with " << "A = " << a << ", B = " << b << ", C = " << c )

    SCAI_REGION( "Storage.CSR.timesMatrix" )

    if ( a.getFormat() != Format::CSR )
    {
        SCAI_UNSUPPORTED( a << ": will be converted to CSR for matrix multiply" )
        auto csrA = convert<CSRStorage<ValueType>>( a );
        csrA.sortRows( false );
        matrixTimesMatrix( alpha, csrA, b, beta, c );
        return;
    }

    if ( b.getFormat() != Format::CSR )
    {
        SCAI_UNSUPPORTED( b << ": will be converted to CSR for matrix multiply" )
        auto csrB = convert<CSRStorage<ValueType>>( b );
        csrB.sortRows( false );
        matrixTimesMatrix( alpha, a, csrB, beta, c );
        return;
    }

    if ( beta != ValueType( 0 ) && c.getFormat() != Format::CSR )
    {
        SCAI_UNSUPPORTED( c << ": CSR temporary required for matrix add" )
        auto csrC = convert<CSRStorage<ValueType>>( c );
        csrC.sortRows( false );
        matrixTimesMatrix( alpha, a, b, beta, csrC );
        return;
    }

    if ( beta == ValueType( 0 ) )
    {
        // this = alpha * a + beta * b 

        matrixTimesMatrixCSR( alpha, static_cast<const CSRStorage<ValueType>&>( a ),
                                     static_cast<const CSRStorage<ValueType>&>( b ) );
    }
    else
    {
        // use a temporary for a * b as this might be aliased to c

        CSRStorage tmp;

        tmp.matrixTimesMatrixCSR( alpha, static_cast<const CSRStorage<ValueType>&>( a ),
                                         static_cast<const CSRStorage<ValueType>&>( b ) );

        matrixPlusMatrixImpl( 1, tmp, beta, static_cast<const CSRStorage<ValueType>&>( c ) );
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void CSRStorage<ValueType>::binaryOp(
    const MatrixStorage<ValueType>& a,
    common::BinaryOp op,
    const MatrixStorage<ValueType>& b )
{
    if ( a.getFormat() != Format::CSR )
    {
        binaryOp( convert<CSRStorage<ValueType>>( a ), op, b );
    }
    else if ( b.getFormat() != Format::CSR )
    {
        binaryOp( a, op, convert<CSRStorage<ValueType>>( b ) );
    }
    else
    {
        binaryOpCSR( static_cast<const CSRStorage<ValueType>&>( a ), op,
                     static_cast<const CSRStorage<ValueType>&>( b ) );
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void CSRStorage<ValueType>::binaryOpCSR(
    const CSRStorage<ValueType>& a,
    common::BinaryOp op,
    const CSRStorage<ValueType>& b )
{
    SCAI_LOG_INFO( logger, "this = a ( " << a.getNumRows() << " x " << a.getNumColumns() << ", nnz = " << a.getNumValues() << " ) " 
                    << op << " b ( " << b.getNumRows() << " x " << b.getNumColumns() << ", nnz = " << b.getNumValues() << " )" )

    if ( &a == this || &b == this )
    {
        // due to alias we use a temporary so that a or b is not damaged during write

        CSRStorage<ValueType> tmp;
        tmp.setContextPtr( getContextPtr() ); 
        tmp.binaryOpCSR( a, op, b );
        swap( tmp ); // safe as tmp will be destroyed afterwards

        return;
    }

    SCAI_ASSERT_EQ_DEBUG( a.getNumRows(), b.getNumRows(), "#rows must be equal for storages in binary op" )
    SCAI_ASSERT_EQ_DEBUG( a.getNumColumns(), b.getNumColumns(), "#cols must be equal for storages in binary op" )

    IndexType m = a.getNumRows();
    IndexType n = a.getNumColumns();

    sparsekernel::CSRUtils::binaryOp( mIA, mJA, mValues, 
                                      a.getIA(), a.getJA(), a.getValues(),
                                      b.getIA(), b.getJA(), b.getValues(),
                                      m, n, op, getContextPtr() );
    
    SCAI_LOG_DEBUG( logger, *this << ": compress by removing zero elements" )

    // by calling the move constructor/assignment we guarantee a consistent storage

    *this = CSRStorage( m, n, std::move( mIA ), std::move( mJA ), std::move( mValues ) );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void CSRStorage<ValueType>::matrixPlusMatrixImpl(
    const ValueType alpha,
    const CSRStorage<ValueType>& a,
    const ValueType beta,
    const CSRStorage<ValueType>& b )
{
    SCAI_LOG_INFO( logger,
                   "this = " << alpha << " * A + " << beta << " * B, with " << "A = " << a << ", B = " << b << ", all are CSR" )

    if ( &a == this || &b == this )
    {
        // due to alias we would get problems with Write/Read access, so use a temporary
        CSRStorage<ValueType> tmp( getContextPtr() );
        tmp.matrixPlusMatrixImpl( alpha, a, beta, b );
        swap( tmp ); // safe as tmp will be destroyed afterwards
        return;
    }

    SCAI_ASSERT_EQ_DEBUG( a.getNumRows(), b.getNumRows(), "serious size mismatch for matrixAdd" )
    SCAI_ASSERT_EQ_DEBUG( a.getNumColumns(), b.getNumColumns(), "serious size mismatch for matrixAdd" )

    IndexType m = a.getNumRows();
    IndexType n = a.getNumColumns();

    sparsekernel::CSRUtils::matrixAdd( mIA, mJA, mValues, 
                                       alpha, a.getIA(), a.getJA(), a.getValues(),
                                       beta, b.getIA(), b.getJA(), b.getValues(),
                                       m, n, getContextPtr() );

    SCAI_LOG_DEBUG( logger, *this << ": compress by removing zero elements" )

    // by calling the move constructor/assignment we guarantee a consistent storage

    *this = CSRStorage( m, n, std::move( mIA ), std::move( mJA ), std::move( mValues ) );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void CSRStorage<ValueType>::matrixTimesMatrixCSR(
    const ValueType alpha,
    const CSRStorage<ValueType>& a,
    const CSRStorage<ValueType>& b )
{
    SCAI_LOG_INFO( logger, *this << ": = " << alpha << " * A * B"
                   << ", with A = " << a << ", B = " << b
                   << ", all are CSR" )

    if ( &a == this || &b == this )
    {
        // due to alias we would get problems with Write/Read access, so use a temporary

        CSRStorage<ValueType> tmp( getContextPtr() );

        tmp.matrixTimesMatrixCSR( alpha, a, b );

        swap( tmp ); // safe as tmp will be destroyed afterwards

        return;
    }

    const IndexType m = a.getNumRows();
    const IndexType n = b.getNumColumns();

    SCAI_ASSERT_EQ_DEBUG( a.getNumColumns(), b.getNumRows(), "serious size mismatch for matrix multiply" )

    const IndexType k = a.getNumColumns();

    // no alias, so reuse arrays of this storage directly

    sparsekernel::CSRUtils::matrixMultiply( mIA, mJA, mValues, alpha,
                                            a.getIA(), a.getJA(), a.getValues(),
                                            b.getIA(), b.getJA(), b.getValues(), 
                                            m, n, k, getContextPtr() );

    *this = CSRStorage( m, n, std::move( mIA ), std::move( mJA ), std::move( mValues ) );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
RealType<ValueType> CSRStorage<ValueType>::l1Norm() const
{
    SCAI_LOG_INFO( logger, *this << ": l1Norm()" )
    return HArrayUtils::l1Norm( mValues, this->getContextPtr() );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
RealType<ValueType> CSRStorage<ValueType>::l2Norm() const
{
    SCAI_LOG_INFO( logger, *this << ": l2Norm()" )
    ContextPtr prefLoc = this->getContextPtr();
    RealType<ValueType> res = HArrayUtils::dotProduct( mValues, mValues, prefLoc );
    res = common::Math::sqrt( res );
    return res;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
RealType<ValueType> CSRStorage<ValueType>::maxNorm() const
{
    // no more checks needed here
    SCAI_LOG_INFO( logger, *this << ": maxNorm()" )

    if ( mValues.size() == 0 )
    {
        return 0;
    }

    static LAMAKernel<UtilKernelTrait::reduce<ValueType> > reduce;
    ContextPtr loc = this->getContextPtr();
    reduce.getSupportedContext( loc );
    ReadAccess<ValueType> csrValues( mValues, loc );
    SCAI_CONTEXT_ACCESS( loc )
    ValueType zero   = 0;
    ValueType maxval = reduce[loc]( csrValues.get(), mValues.size(), zero, BinaryOp::ABS_MAX );
    return maxval;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
RealType<ValueType> CSRStorage<ValueType>::maxDiffNorm( const MatrixStorage<ValueType>& other ) const
{
    SCAI_REGION( "Storage.CSR.maxDiffNorm" )

    SCAI_ASSERT_EQ_ERROR( getNumRows(), other.getNumRows(), "row size mismatch for maxDiffNorm" )
    SCAI_ASSERT_EQ_ERROR( getNumColumns(), other.getNumColumns(), "col size mismatch for maxDiffNorm" )

    SCAI_LOG_INFO( logger, *this << ": maxDiffNorm( " << other << " )" )

    if ( other.getFormat() == Format::CSR )
    {
        return maxDiffNormImpl( static_cast<const CSRStorage<ValueType>&>( other ) );
    }
    else
    {
        SCAI_UNSUPPORTED( other << ": converted to " << typeName() << " for maxDiffNorm" )

        return maxDiffNormImpl( convert<CSRStorage<ValueType>>( other ) );
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
RealType<ValueType> CSRStorage<ValueType>::maxDiffNormImpl( const CSRStorage<ValueType>& other ) const
{
    // no more checks needed here
    SCAI_LOG_INFO( logger, *this << ": maxDiffNormImpl( " << other << " )" )

    if ( getNumRows() == 0 )
    {
        return static_cast<ValueType>( 0.0 );
    }

    bool sorted = mSortedRows && other.mSortedRows && ( mDiagonalProperty == other.mDiagonalProperty );
    static LAMAKernel<CSRKernelTrait::absMaxDiffVal<ValueType> > absMaxDiffVal;
    ContextPtr loc = this->getContextPtr();
    absMaxDiffVal.getSupportedContext( loc );
    ReadAccess<IndexType> csrIA1( mIA, loc );
    ReadAccess<IndexType> csrJA1( mJA, loc );
    ReadAccess<ValueType> csrValues1( mValues, loc );
    ReadAccess<IndexType> csrIA2( other.mIA, loc );
    ReadAccess<IndexType> csrJA2( other.mJA, loc );
    ReadAccess<ValueType> csrValues2( other.mValues, loc );
    SCAI_CONTEXT_ACCESS( loc )
    ValueType maxval = absMaxDiffVal[loc] ( getNumRows(), sorted, csrIA1.get(), csrJA1.get(), csrValues1.get(), csrIA2.get(),
                                            csrJA2.get(), csrValues2.get() );
    return maxval;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
CSRStorage<ValueType>* CSRStorage<ValueType>::copy() const
{
    return new CSRStorage<ValueType>( *this );
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
void CSRStorage<ValueType>::buildSparseRowSizes( HArray<IndexType>& rowSizes ) const
{
    SCAI_LOG_DEBUG( logger, "copy nnz for each row in HArray" );
    WriteOnlyAccess<IndexType> writeRowSizes( rowSizes, getNumRows() );
    ReadAccess<IndexType> csrIA( mIA );
    OpenMPCSRUtils::offsets2sizes( writeRowSizes.get(), csrIA.get(), getNumRows() );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void CSRStorage<ValueType>::buildSparseRowData(
    HArray<IndexType>& sparseJA,
    HArray<ValueType>& sparseValues ) const
{
    SCAI_LOG_INFO( logger, *this << ": build sparse row data" );
    // for CSR format we can just copy arrays with column indexes and data values
    sparseJA = mJA;
    sparseValues = mValues;
}

/* ========================================================================= */

template<typename ValueType>
void CSRStorage<ValueType>::fillCOO(
    hmemo::HArray<IndexType> cooIA,
    hmemo::HArray<IndexType> cooJA,
    hmemo::HArray<ValueType> cooValues,
    const common::BinaryOp op )
{
    using sparsekernel::COOUtils;

    SCAI_ASSERT_EQ_ERROR( cooIA.size(), cooJA.size(), "COO data: cooIA and cooJA must have same size" )
    SCAI_ASSERT_EQ_ERROR( cooIA.size(), cooValues.size(), "COO data: cooIA and cooValues must have same size" )

    // convert the COO data to CSR, so it is sorted and filling can be done in parallel

    COOUtils::normalize( cooIA, cooJA, cooValues, op, getContextPtr() );

    HArray<IndexType> csrIA;

    COOUtils::convertCOO2CSR( csrIA, cooIA, getNumRows(), getContextPtr() );

    CSRStorage<ValueType> csr1( getNumRows(), getNumColumns(), std::move( csrIA ), std::move( cooJA ), std::move( cooValues ) );

    if ( op == common::BinaryOp::COPY )
    {
        csr1.writeToFile( "COO_copy.mtx" );
    }
    else 
    {
        csr1.writeToFile( "COO_add.mtx" );
    }

    sortRows( false );

    this->writeToFile( "CSR.mtx" );

    binaryOpCSR( *this, op, csr1 );

    if ( op == common::BinaryOp::COPY )
    {
        this->writeToFile( "CSR_copy.mtx" );
    }
    else 
    {
        this->writeToFile( "CSR_add.mtx" );
    }
}

/* ========================================================================= */
/*  Static fatory methods and related virtual methods                        */
/* ========================================================================= */

template<typename ValueType>
std::string CSRStorage<ValueType>::initTypeName()
{
    std::stringstream s;
    s << std::string( "CSRStorage<" ) << common::getScalarType<ValueType>() << std::string( ">" );
    return s.str();
}

template<typename ValueType>
const char* CSRStorage<ValueType>::typeName()
{
    static const std::string s = initTypeName();
    return  s.c_str();
}

template<typename ValueType>
const char* CSRStorage<ValueType>::getTypeName() const
{
    return typeName();
}

template<typename ValueType>
MatrixStorageCreateKeyType CSRStorage<ValueType>::createValue()
{
    return MatrixStorageCreateKeyType( Format::CSR, common::getScalarType<ValueType>() );
}

template<typename ValueType>
MatrixStorageCreateKeyType CSRStorage<ValueType>::getCreateValue() const
{
    return createValue();
}

template<typename ValueType>
_MatrixStorage* CSRStorage<ValueType>::create()
{
    return new CSRStorage<ValueType>();
}

/* ========================================================================= */
/*       Template Instantiations                                             */
/* ========================================================================= */

SCAI_COMMON_INST_CLASS( CSRStorage, SCAI_NUMERIC_TYPES_HOST )

#define CSR_STORAGE_INST_LVL2( ValueType, OtherValueType )    \
    template void CSRStorage<ValueType>::setCSRDataImpl(      \
            const IndexType,                                  \
            const IndexType,                                  \
            const hmemo::HArray<IndexType>&,                  \
            const hmemo::HArray<IndexType>&,                  \
            const hmemo::HArray<OtherValueType>&,             \
            const hmemo::ContextPtr );                        \
    template void CSRStorage<ValueType>::buildCSR(            \
            hmemo::HArray<IndexType>&,                        \
            hmemo::HArray<IndexType>*,                        \
            hmemo::HArray<OtherValueType>*,                   \
            const hmemo::ContextPtr ) const;                  

#define CSR_STORAGE_INST_LVL1( ValueType )                                                                                  \
    SCAI_COMMON_LOOP_LVL2( ValueType, CSR_STORAGE_INST_LVL2, SCAI_NUMERIC_TYPES_HOST )

SCAI_COMMON_LOOP( CSR_STORAGE_INST_LVL1, SCAI_NUMERIC_TYPES_HOST )

#undef CSR_STORAGE_INST_LVL2
#undef CSR_STORAGE_INST_LVL1


} /* end namespace lama */

} /* end namespace scai */
