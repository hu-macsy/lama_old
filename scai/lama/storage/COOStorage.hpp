/**
 * @file COOStorage.hpp
 *
 * @license
 * Copyright (c) 2009-2015
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
 * @brief Definition of a matrix storage class for a (non-distributed) sparse matrix
 *        represented by coordinate lists.
 * @author Paul Börding, Thomas Brandes
 * @date 11.06.2011
 */
#pragma once

// for dll_import
#include <scai/common/config.hpp>

// base classes
#include <scai/lama/storage/CRTPMatrixStorage.hpp>

// internal scai libraries
#include <scai/sparsekernel/CSRKernelTrait.hpp>
#include <scai/sparsekernel/COOKernelTrait.hpp>
#include <scai/utilskernel/LAMAKernel.hpp>
#include <scai/utilskernel/UtilKernelTrait.hpp>
#include <scai/common/unique_ptr.hpp>

namespace scai
{

namespace lama
{

/** Matrix storage format for a COO sparse matrix.
 *
 * COO stores a list of (row, column, value) tuples. For efficiency reasons,
 * three separate arrays are used instead of a single array of triples.
 *
 * The entries might be sorted (e.g. by row index or by column index) to improve random access times.
 * The current strategy is to sort them row-wise and then column-wise so that access patterns
 * are the same as in the CSR and ELL format. But the actual order is not important for the
 * correctness of the implemented methods.
 *
 * If the diagonal property flag is set, the first numRows entries of the COO data will
 * contain the diagonal elements.
 *
 * @tparam ValueType is the value type of the matrix values.
 *
 * Note: Copy constructor and operator= are supported by default implementations,
 *       copies of HArray will be deep copies.
 */

template<typename ValueType>
class COMMON_DLL_IMPORTEXPORT COOStorage:
    public CRTPMatrixStorage<COOStorage<ValueType>,ValueType>,
    public _MatrixStorage::Register<COOStorage<ValueType> >    // register at factory
{
public:

    typedef ValueType StorageValueType;

    /** get typename of the matrix storage format. */

    static const char* typeName();

    /**
     * @brief Creates a sparse matrix with all values set to zero.
     *
     * @param[in] numRows       the number of rows of the matrix
     * @param[in] numColumns    the number of columns of the matrix
     */
    COOStorage( const IndexType numRows, const IndexType numColumns );

    /**
     * @brief Creates a sparse matrix with the coo array
     */
    COOStorage(
        const IndexType numRows,
        const IndexType numColumns,
        const hmemo::HArray<IndexType>& ia,
        const hmemo::HArray<IndexType>& ja,
        const hmemo::HArray<ValueType>& values );

    /** Default constructor, same as COOStorage(0, 0). */

    COOStorage();

    /** Default copy constructor is overridden */

    COOStorage( const COOStorage<ValueType>& other );

    /** Copy constructor that handles also type and format conversion. */

    explicit COOStorage( const _MatrixStorage& other )
    {
        assign( other );
    }

    /** Copy constructor can take any matrix storage or context. */

    COOStorage( const _MatrixStorage& other, const hmemo::ContextPtr context )
    {
        setContextPtr( context );
        assign( other );
    }

    /** Default assignment operator is overridden */

    COOStorage<ValueType>& operator=( const COOStorage<ValueType>& other )
    {
        assign( other );
        return *this;
    }

    COOStorage<ValueType>& operator=( const _MatrixStorage& other )
    {
        assign( other );
        return *this;
    }

    virtual void clear();

    /** Destructor of COO sparse matrix. */

    virtual ~COOStorage();

    /** Test the storage data for inconsistencies.
     *
     *  @throw Exception in case of any inconsistency.
     */
    void check( const char* msg ) const;

    /** Getter routine for the enum value that stands for this format. */

    virtual Format::MatrixStorageFormat getFormat() const;

    /** Resize of a zero matrix.
     *
     *  \code
     *     COOStorage a(0, 0);
     *     a.allocate(100, 100)
     *  \endcode
     */

    void allocate( const IndexType numRows, const IndexType numColumns );

    /** Implementation of pure method of class MatrixStorage. */

    virtual void purge();

    /** Implemenation of pure method of class MatrixStorage. */

    virtual void print( std::ostream& ) const;

    /** Implementation of pure method for COO matrix. */

    virtual void setIdentity( const IndexType size );

    /** Help routine to build any kind of CSR storage. */

    template<typename OtherValueType>
    void buildCSR(
        hmemo::HArray<IndexType>& ia,
        hmemo::HArray<IndexType>* ja,
        hmemo::HArray<OtherValueType>* values,
        const hmemo::ContextPtr context ) const;

    /**
     * @brief fills COO sparse matrix by coo sparse data.
     *
     * @param[in] numRows    number of rows
     * @param[in] numColumns number of columns
     * @param[in] numValues  the number of stored elements in the matrix
     * @param[in] ia         row indexes for each value
     * @param[in] ja         column indexes for each value
     * @param[in] values     data values of the corresponding matrix elements
     */
    void setCOOData(
        const IndexType numRows,
        const IndexType numColumns,
        const IndexType numValues,
        const hmemo::HArray<IndexType>& ia,
        const hmemo::HArray<IndexType>& ja,
        const hmemo::_HArray& values );

    /**
     * @brief fills COO sparse matrix by csr sparse data.
     *
     * @param[in] numRows    number of rows
     * @param[in] numColumns number of columns
     * @param[in] numValues  the number of stored elements in the matrix
     * @param[in] ia         row pointer of the input csr sparse matrix
     * @param[in] ja         column indexes of the input csr sparse matrix
     * @param[in] values     the data values of the input csr sparse matrix
     * @param[in] loc        is the context where filling takes place
     */
    template<typename OtherValueType>
    void setCSRDataImpl(
        const IndexType numRows,
        const IndexType numColumns,
        const IndexType numValues,
        const hmemo::HArray<IndexType>& ia,
        const hmemo::HArray<IndexType>& ja,
        const hmemo::HArray<OtherValueType>& values,
        const hmemo::ContextPtr loc ) __attribute__( ( noinline ) );

    /* Print relevant information about matrix storage format. */

    virtual void writeAt( std::ostream& stream ) const;

    /** Getter for member variable IA, only const reference. */

    const hmemo::HArray<IndexType>& getIA() const;

    /** Getter for member variable JA, only const reference. */

    const hmemo::HArray<IndexType>& getJA() const;

    /** Getter for member variable values, only const reference. */

    const hmemo::HArray<ValueType>& getValues() const;

    /** Getter for the number of stored values. */

    virtual IndexType getNumValues() const;

    /** Template method for getting row. */

    template<typename OtherType>
    void getRowImpl( hmemo::HArray<OtherType>& row, const IndexType i ) const __attribute( ( noinline ) );

    /** This method returns the diagonal
     *
     * @param[in] diagonal  is the destination array
     *
     * Calculations are dependent to the diagonal property
     */
    template<typename OtherType>
    void getDiagonalImpl( hmemo::HArray<OtherType>& diagonal ) const __attribute( ( noinline ) );

    /** Template version used for virtual routine setDiagonal with known value type. */

    template<typename OtherType>
    void setDiagonalImpl( const hmemo::HArray<OtherType>& diagonal ) __attribute( ( noinline ) );

    /** Implementation of pure method. */

    virtual void setDiagonalImpl( const ValueType value );

    /** Implementation of MatrixStorage::copy for derived class. */

    virtual COOStorage* copy() const;

    /** Implementation of MatrixStorage::newMatrixStorage for derived class. */

    virtual COOStorage* newMatrixStorage() const;

    /******************************************************************
     *  Matrix times Vector                                            *
     ******************************************************************/

    /** Implementation of MatrixStorage::matrixTimesVector for COO */

    virtual void matrixTimesVector(
        hmemo::HArray<ValueType>& result,
        const ValueType alpha,
        const hmemo::HArray<ValueType>& x,
        const ValueType beta,
        const hmemo::HArray<ValueType>& y ) const;

    /** Implementation of MatrixStorage::vectorTimesMatrix for COO */
    /** since 1.0.1 */

    virtual void vectorTimesMatrix(
        hmemo::HArray<ValueType>& result,
        const ValueType alpha,
        const hmemo::HArray<ValueType>& x,
        const ValueType beta,
        const hmemo::HArray<ValueType>& y ) const;

    /** Implementation of MatrixStorage::matrixTimesVectorAsync for COO */

    virtual tasking::SyncToken* matrixTimesVectorAsync(
        hmemo::HArray<ValueType>& result,
        const ValueType alpha,
        const hmemo::HArray<ValueType>& x,
        const ValueType beta,
        const hmemo::HArray<ValueType>& y ) const;

    /** Implementation of MatrixStorage::vectorTimesMatrixAsync for CSR */
    /** since 1.0.1 */

    virtual tasking::SyncToken* vectorTimesMatrixAsync(
        hmemo::HArray<ValueType>& result,
        const ValueType alpha,
        const hmemo::HArray<ValueType>& x,
        const ValueType beta,
        const hmemo::HArray<ValueType>& y ) const;

    /** solution = xxx */

    virtual void jacobiIterate(
        hmemo::HArray<ValueType>& solution,
        const hmemo::HArray<ValueType>& oldSolution,
        const hmemo::HArray<ValueType>& rhs,
        const ValueType omega ) const;

    /******************************************************************
     *  Scaling of elements in a matrix                                *
     ******************************************************************/

    /** Template version used for virtual routine scale with known value type. */

    template<typename OtherType>
    void scaleImpl( const hmemo::HArray<OtherType>& values ) __attribute( ( noinline ) );

    /** Implementation of pure method.  */

    void scaleImpl( const ValueType value );

    /** Implementation of pure method.  */

    void conj();

    /** Implementation for MatrixStorage::l1Norm */

    virtual ValueType l1Norm() const;

    /** Implementation for MatrixStorage::l2Norm */

    virtual ValueType l2Norm() const;

    /** Implementation for MatrixStorage::maxNorm */

    virtual ValueType maxNorm() const;

    /** Get a value of the matrix.
     *
     * @param[in] i is the row index, 0 <= i < mNumRows
     * @param[in] j is the colum index, 0 <= j < mNumRows
     *
     * Out-of-range check is enabled for DEBUG version.
     */

    ValueType getValue( const IndexType i, const IndexType j ) const;

    /** Initiate an asynchronous data transfer to a specified location. */

    void prefetch( const hmemo::ContextPtr location ) const;

    /** Will wait for all outstanding asynchronous data transfers. */

    void wait() const;

    /** Swaps this with other.
     * @param[in,out] other the COOStorage to swap this with
     */
    void swap( COOStorage<ValueType>& other );

    virtual size_t getMemoryUsageImpl() const;

    using MatrixStorage<ValueType>::assign;
    using MatrixStorage<ValueType>::prefetch;
    using MatrixStorage<ValueType>::getContextPtr;
    using MatrixStorage<ValueType>::setContextPtr;

protected:

    using MatrixStorage<ValueType>::mNumRows;
    using MatrixStorage<ValueType>::mNumColumns;
    using MatrixStorage<ValueType>::mDiagonalProperty;
    using MatrixStorage<ValueType>::mRowIndexes;
    using MatrixStorage<ValueType>::mCompressThreshold;

    IndexType mNumValues; //!< number of non-zero values (+ optionally zeros in diagonal)

    hmemo::HArray<IndexType> mIA; //!< row indices, size is mNumValues
    hmemo::HArray<IndexType> mJA; //!< column indices, size is mNumValues
    hmemo::HArray<ValueType> mValues; //!< non-zero values (+ optionally zeros in diagonal), size is mNumValues

private:

    /** result += alpha * (*this) * x */

    tasking::SyncToken* incGEMV( hmemo::HArray<ValueType>& result,
                                 const ValueType alpha,
                                 const hmemo::HArray<ValueType>& x,
                                 bool async ) const;

    /** result += alpha * x * (*this) */

    tasking::SyncToken* incGEVM( hmemo::HArray<ValueType>& result,
                                 const ValueType alpha,
                                 const hmemo::HArray<ValueType>& x,
                                 bool async ) const;

    /** Function that checks the diagonal property of the COO matrix. If on
     *  the first column index for each row is the row (diagonal element).
     */
    virtual bool checkDiagonalProperty() const;

    // values are stored in row-major order

    // TODO: index at irow / icolumn, not supported yet

    // inline IndexType index(IndexType irow, IndexType icolumn) const { return icolumn * mNumRows + irow; }

    SCAI_LOG_DECL_STATIC_LOGGER( logger )

    static std::string initTypeName();

public:

    // static create method that will be used to register at MatrixStorage factory

    static _MatrixStorage* create();

    // key for factory

    static MatrixStorageCreateKeyType createValue();
};

template<typename ValueType>
template<typename OtherValueType>
void COOStorage<ValueType>::buildCSR(
    hmemo::HArray<IndexType>& ia,
    hmemo::HArray<IndexType>* ja,
    hmemo::HArray<OtherValueType>* values,
    const hmemo::ContextPtr preferredLoc ) const
{
    // multiple kernel routines needed

    static utilskernel::LAMAKernel<sparsekernel::CSRKernelTrait::sizes2offsets> sizes2offsets;
    static utilskernel::LAMAKernel<sparsekernel::COOKernelTrait::getCSRSizes> getCSRSizes;
    static utilskernel::LAMAKernel<sparsekernel::COOKernelTrait::getCSRValues<ValueType, OtherValueType> > getCSRValues;

    // do it where all routines are avaialble

    hmemo::ContextPtr loc = preferredLoc;

    sizes2offsets.getSupportedContext( loc, getCSRSizes, getCSRValues );

    SCAI_CONTEXT_ACCESS( loc )

    hmemo::WriteOnlyAccess<IndexType> csrIA( ia, loc, mNumRows + 1 );
    hmemo::ReadAccess<IndexType> cooIA( mIA, loc );

    getCSRSizes[loc]( csrIA.get(), mNumRows, mNumValues, cooIA.get() );

    if( ja == NULL || values == NULL )
    {
        csrIA.resize( mNumRows );
        return;
    }

    IndexType numValues = sizes2offsets[loc]( csrIA.get(), mNumRows );

    SCAI_ASSERT_EQUAL_DEBUG( mNumValues, numValues )

    hmemo::ReadAccess<IndexType> cooJA( mJA, loc );
    hmemo::ReadAccess<ValueType> cooValues( mValues, loc );

    hmemo::WriteOnlyAccess<IndexType> csrJA( *ja, loc, numValues );
    hmemo::WriteOnlyAccess<OtherValueType> csrValues( *values, loc, numValues );

    getCSRValues[loc]( csrJA.get(), csrValues.get(), csrIA.get(),
                       mNumRows, mNumValues, cooIA.get(), cooJA.get(), cooValues.get() );
}


template<typename ValueType>
template<typename OtherValueType>
void COOStorage<ValueType>::setCSRDataImpl(
    const IndexType numRows,
    const IndexType numColumns,
    const IndexType numValues,
    const hmemo::HArray<IndexType>& ia,
    const hmemo::HArray<IndexType>& ja,
    const hmemo::HArray<OtherValueType>& values,
    const hmemo::ContextPtr prefLoc )
{
    SCAI_LOG_DEBUG( logger, "set CSR data " << numRows << " x " << numColumns << ", nnz = " << numValues )

    if ( ia.size() == numRows )
    {
        // offset array required

        hmemo::HArray<IndexType> offsets;

        IndexType total = _MatrixStorage::sizes2offsets( offsets, ia, prefLoc );

        SCAI_ASSERT_EQUAL( numValues, total, "sizes do not sum to number of values" );

        setCSRDataImpl( numRows, numColumns, numValues, offsets, ja, values, prefLoc );

        return;
    }

    SCAI_ASSERT_EQUAL_DEBUG( numRows + 1, ia.size() )
    SCAI_ASSERT_EQUAL_DEBUG( numValues, ja.size() )
    SCAI_ASSERT_EQUAL_DEBUG( numValues, values.size() )

    hmemo::ContextPtr loc = prefLoc;

    // ReadAccess<IndexType> csrJA( ja, loc );
    // ReadAccess<OtherValueType> csrValues( values, loc );

    mNumRows = numRows;
    mNumColumns = numColumns;

    // check if input csr data has the diagonal property and inherit it

    int numDiagonals = std::min( numRows, numColumns );

    {
        SCAI_LOG_DEBUG( logger,
                        "check CSR data " << numRows << " x " << numColumns << ", nnz = " << numValues << " for diagonal property, #diagonals = " << numDiagonals )

        static utilskernel::LAMAKernel<sparsekernel::CSRKernelTrait::hasDiagonalProperty> hasDiagonalProperty;

        hmemo::ContextPtr loc = this->getContextPtr();

        hasDiagonalProperty.getSupportedContext( loc );

        hmemo::ReadAccess<IndexType> csrIA( ia, loc );
        hmemo::ReadAccess<IndexType> csrJA( ja, loc );

        SCAI_CONTEXT_ACCESS( loc )

        mDiagonalProperty = hasDiagonalProperty[loc] ( numDiagonals, csrIA.get(), csrJA.get() );
    }

    if( !mDiagonalProperty )
    {
        numDiagonals = 0; // do not store diagonal data at the beginning in COO data
    }

    mNumValues = numValues;

    SCAI_LOG_DEBUG( logger,
                    "input csr data with " << mNumValues << "entries,  has diagonal property = " << mDiagonalProperty )

    {
        static utilskernel::LAMAKernel<sparsekernel::COOKernelTrait::offsets2ia> offsets2ia;

        hmemo::ContextPtr loc = this->getContextPtr();
        offsets2ia.getSupportedContext( loc );

        hmemo::ReadAccess<IndexType> csrIA( ia, loc );
        hmemo::WriteOnlyAccess<IndexType> cooIA( mIA, loc, mNumValues );

        SCAI_CONTEXT_ACCESS( loc )

        offsets2ia[loc]( cooIA.get(), mNumValues, csrIA.get(), mNumRows, numDiagonals );
    }

    {
        static utilskernel::LAMAKernel<sparsekernel::COOKernelTrait::setCSRData<IndexType, IndexType> > setCSRData;

        hmemo::ContextPtr loc = this->getContextPtr();   // preferred location
        setCSRData.getSupportedContext( loc );    // supported location

        hmemo::ReadAccess<IndexType> csrIA( ia, loc );
        hmemo::ReadAccess<IndexType> csrJA( ja, loc );
        hmemo::WriteOnlyAccess<IndexType> cooJA( mJA, loc, mNumValues );

        SCAI_CONTEXT_ACCESS( loc )

        setCSRData[loc]( cooJA.get(), csrJA.get(), numValues, csrIA.get(), mNumRows, numDiagonals );
    }

    {
        static utilskernel::LAMAKernel<sparsekernel::COOKernelTrait::setCSRData<ValueType, OtherValueType> > setCSRData;

        hmemo::ContextPtr loc = this->getContextPtr();   // preferred location
        setCSRData.getSupportedContext( loc );    // supported location

        hmemo::ReadAccess<IndexType> csrIA( ia, loc );
        hmemo::ReadAccess<OtherValueType> csrValues( values, loc );
        hmemo::WriteOnlyAccess<ValueType> cooValues( mValues, loc, mNumValues );

        SCAI_CONTEXT_ACCESS( loc )

        setCSRData[loc]( cooValues.get(), csrValues.get(), numValues, csrIA.get(), mNumRows, numDiagonals );
    }
}

template<typename ValueType>
template<typename OtherType>
void COOStorage<ValueType>::setDiagonalImpl( const hmemo::HArray<OtherType>& diagonal )
{
    const IndexType numDiagonalElements = std::min( mNumColumns, mNumRows );

    static utilskernel::LAMAKernel<utilskernel::UtilKernelTrait::set<ValueType, OtherType> > set;

    hmemo::ContextPtr loc = this->getContextPtr();
    set.getSupportedContext( loc );

    hmemo::ReadAccess<OtherType> rDiagonal( diagonal, loc );
    hmemo::WriteAccess<ValueType> wValues( mValues, loc );

    SCAI_CONTEXT_ACCESS( loc )

    // diagonal elements are the first entries of mValues

    set[loc]( wValues.get(), rDiagonal.get(), numDiagonalElements, common::reduction::COPY );
}

template<typename ValueType>
template<typename OtherType>
void COOStorage<ValueType>::getDiagonalImpl( hmemo::HArray<OtherType>& diagonal ) const
{
    // diagional[0:numDiagonalElements] = mValues[0:numDiagonalElements]
    // Note: using HArrayUtils::setArray not possible, as we only need part of mValues

    const IndexType numDiagonalElements = std::min( mNumColumns, mNumRows );

    static utilskernel::LAMAKernel<utilskernel::UtilKernelTrait::set<OtherType, ValueType> > set;

    hmemo::ContextPtr loc = this->getContextPtr();
    set.getSupportedContext( loc );

    hmemo::WriteOnlyAccess<OtherType> wDiagonal( diagonal, loc, numDiagonalElements );
    hmemo::ReadAccess<ValueType> rValues( mValues, loc );

    SCAI_CONTEXT_ACCESS( loc )

    // diagonal elements are the first entries of mValues

    set[loc]( wDiagonal.get(), rValues.get(), numDiagonalElements, common::reduction::COPY );
}

template<typename ValueType>
template<typename OtherType>
void COOStorage<ValueType>::getRowImpl( hmemo::HArray<OtherType>& row, const IndexType i ) const
{
    SCAI_ASSERT_DEBUG( i >= 0 && i < mNumRows, "row index " << i << " out of range" )

                hmemo::ContextPtr hostContext = hmemo::Context::getHostPtr();

    hmemo::WriteOnlyAccess<OtherType> wRow( row, mNumColumns );

    const hmemo::ReadAccess<IndexType> ia( mIA, hostContext );
    const hmemo::ReadAccess<IndexType> ja( mJA, hostContext );
    const hmemo::ReadAccess<ValueType> values( mValues, hostContext );

    // ToDo: OpenMP parallelization, interface

    for( IndexType j = 0; j < mNumColumns; ++j )
    {
        wRow[j] = static_cast<OtherType>(0.0);
    }

    for( IndexType kk = 0; kk < mNumValues; ++kk )
    {
        if( ia[kk] != i )
        {
            continue;
        }

        wRow[ja[kk]] = static_cast<OtherType>( values[kk] );
    }
}

template<typename ValueType>
void COOStorage<ValueType>::scaleImpl( const ValueType value )
{
    // multiply value with each entry of mValues

    utilskernel::HArrayUtils::setScalar( mValues, value, common::reduction::MULT, this->getContextPtr() );
}

} /* end namespace lama */

} /* end namespace scai */
