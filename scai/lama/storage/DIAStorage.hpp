/**
 * @file DIAStorage.hpp
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
 * @brief Definition of a structure for a (non-distributed) DIA sparse matrix.
 * @author Thomas Brandes
 * @date 11.06.2011
 */

#pragma once

// for dll_import
#include <scai/common/config.hpp>

// base classes
#include <scai/lama/storage/CRTPMatrixStorage.hpp>

namespace scai
{

namespace lama
{

/** @brief Storage format for a DIA (compressed storage row) sparse matrix.
 *
 *  The DIA format has the following data:
 *
 *  - number of diagonals
 *  - offsets of the diagonals
 *  - values of the diagonals
 *
 *  @tparam ValueType is the value type of the matrix values.
 */
template<typename ValueType>
class COMMON_DLL_IMPORTEXPORT DIAStorage:
    public CRTPMatrixStorage<DIAStorage<ValueType>,ValueType>,
    public _MatrixStorage::Register<DIAStorage<ValueType> >    // register at factory
{
public:

    typedef ValueType StorageValueType;

    /** Getter of type name of the matrix storage format. */

    static const char* typeName();

    /**
     * @brief Creates a sparse matrix with all values set to zero.
     *
     * @param[in] numRows       the number of rows of the matrix
     * @param[in] numColumns    the number of columns of the matrix
     */
    DIAStorage( const IndexType numRows, const IndexType numColumns );

    /** Default constructor, same as DIAStorage(0, 0). */

    DIAStorage();

    /** Constructor for DIA storage by corresponding arrays. */

    DIAStorage(
        const IndexType numRows,
        const IndexType numColumns,
        const IndexType numDiagonals,
        const hmemo::HArray<IndexType>& offsets,
        const hmemo::HArray<ValueType>& values );

    /** Default copy constructor is overridden */

    DIAStorage( const DIAStorage<ValueType>& other );

    /** Copy constructor that handles also type and format conversion. */

    explicit DIAStorage( const _MatrixStorage& other )
    {
        assign( other );
    }

    /** Copy constructor can take any matrix storage or context. */

    DIAStorage( const _MatrixStorage& other, const hmemo::ContextPtr context )
    {
        setContextPtr( context );
        assign( other );
    }

    /** Default assignment operator is overridden */

    DIAStorage<ValueType>& operator=( const DIAStorage<ValueType>& other )
    {
        assign( other );
        return *this;
    }

    DIAStorage<ValueType>& operator=( const _MatrixStorage& other )
    {
        assign( other );
        return *this;
    }

    /** Implementation of MatrixStorage::copy for derived class. */

    virtual DIAStorage* copy() const;

    /** Implementation of _MatrixStorage::clear  */

    virtual void clear();

    /** Destructor of DIA sparse matrix. */

    virtual ~DIAStorage();

    /** Test the storage data for inconsistencies.
     *
     *  @throw Exception in case of any inconsistency.
     */
    void check( const char* msg ) const;

    /** Getter routine for the enum value that stands for this format. */

    virtual MatrixStorageFormat getFormat() const;

    /** Resize of a zero matrix.
     *
     *  \code
     *     DIAStorage a(0, 0);
     *     a.allocate(100, 100)
     *  \endcode
     */

    void allocate( const IndexType numRows, const IndexType numColumns );

    /** Implementation of pure method of class MatrixStorage. */

    virtual void purge();

    /**  Implementation of pure method for DIA storage format.  */

    virtual void setIdentity( const IndexType size );

    /** General routine to build any kind of CSR storage.
     *
     *  @param[out] ia is the CSR offset array
     *  @param[out] ja is the array with the column indexes (optional)
     *  @param[out] values is the array with the non-zero matrix values (optional)
     *  @param[in]  loc is the Context where conversion should be done
     */
    template<typename OtherValueType>
    void buildCSR(
        hmemo::HArray<IndexType>& ia,
        hmemo::HArray<IndexType>* ja,
        hmemo::HArray<OtherValueType>* values,
        const hmemo::ContextPtr loc ) const;

    /**
     * @brief fills DIA matrix storage by csr sparse data.
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
        const hmemo::ContextPtr loc );

    /** Implementation of MatrixStorage::matrixTimesVector for DIA */

    virtual void matrixTimesVector(
        hmemo::HArray<ValueType>& result,
        const ValueType alpha,
        const hmemo::HArray<ValueType>& x,
        const ValueType beta,
        const hmemo::HArray<ValueType>& y ) const;

    /** Implementation of MatrixStorage::vectorTimesMatrix for DIA */

    virtual void vectorTimesMatrix(
        hmemo::HArray<ValueType>& result,
        const ValueType alpha,
        const hmemo::HArray<ValueType>& x,
        const ValueType beta,
        const hmemo::HArray<ValueType>& y ) const;

    /** Implementation of MatrixStorage::matrixTimesVectorAsync for DIA */

    virtual tasking::SyncToken* matrixTimesVectorAsync(
        hmemo::HArray<ValueType>& result,
        const ValueType alpha,
        const hmemo::HArray<ValueType>& x,
        const ValueType beta,
        const hmemo::HArray<ValueType>& y ) const;

    /** Implementation of MatrixStorage::vectorTimesMatrixAsync for DIA */

    virtual tasking::SyncToken* vectorTimesMatrixAsync(
        hmemo::HArray<ValueType>& result,
        const ValueType alpha,
        const hmemo::HArray<ValueType>& x,
        const ValueType beta,
        const hmemo::HArray<ValueType>& y ) const;

    /** Implementation of MatrixStorage::jacobiIterate for DIA */

    virtual void jacobiIterate(
        hmemo::HArray<ValueType>& solution,
        const hmemo::HArray<ValueType>& oldSolution,
        const hmemo::HArray<ValueType>& rhs,
        const ValueType omega ) const;

    /* Print relevant information about matrix storage format. */

    virtual void writeAt( std::ostream& stream ) const;

    /** Getter routine for the matrix offsets (read-only). */

    const hmemo::HArray<IndexType>& getOffsets() const;

    /** Getter routine for the matrix values (read-only). */

    const hmemo::HArray<ValueType>& getValues() const;

    /** Getter routine for the number of diagonals. */

    IndexType getNumDiagonals() const;

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

    /** This method replaces the diagonal by an array of new values
     *
     * @param[in] diagonal  is the array with the new values
     *
     * Calculations are dependent to the diagonal property
     */
    template<typename OtherType>
    void setDiagonalImpl( const hmemo::HArray<OtherType>& diagonal ) __attribute( ( noinline ) );

    /** This method replaces the diagonal by a diagonal value
     *
     * @param[in] value  is the source value
     *
     * Calculations are dependent to the diagonal property
     */
    void setDiagonalImpl( const ValueType value );

    /******************************************************************
     *  Scaling of elements in a matrix                                *
     ******************************************************************/

    /** Template version used for virtual routine scale with known value type. */

    template<typename OtherType>
    void scaleImpl( const hmemo::HArray<OtherType>& values ) __attribute( ( noinline ) );

    /** Implementation of pure method.  */

    virtual void scaleImpl( const ValueType value );

    /** Implementation of pure method.  */

    virtual void conj();

    /** Implementation for MatrixStorage::l1Norm */

    virtual ValueType l1Norm() const;

    /** Implementation for MatrixStorage::l2Norm */

    virtual ValueType l2Norm() const;

    /** Implementation for MatrixStorage::maxNorm */

    virtual ValueType maxNorm() const;

    /** Implementation of pure method. */

    ValueType getValue( const IndexType i, const IndexType j ) const;

    /** Initiate an asynchronous data transfer to a specified location. */

    void prefetch( const hmemo::ContextPtr location ) const;

    /** Will wait for all outstanding asynchronous data transfers. */

    void wait() const;

    /** Swaps this with other.
     * @param[in,out] other the DIAStorage to swap this with
     */
    void swap( DIAStorage<ValueType>& other );

    virtual size_t getMemoryUsageImpl() const;

    /** print matrix on cout, helpful for debug. */

    void print() const;

    using MatrixStorage<ValueType>::prefetch;
    using MatrixStorage<ValueType>::getContextPtr;
    using MatrixStorage<ValueType>::assign;
    using MatrixStorage<ValueType>::setContextPtr;

protected:

    using MatrixStorage<ValueType>::mNumRows;
    using MatrixStorage<ValueType>::mNumColumns;
    using MatrixStorage<ValueType>::mDiagonalProperty;
    using MatrixStorage<ValueType>::mRowIndexes;
    using MatrixStorage<ValueType>::mCompressThreshold;

    IndexType mNumDiagonals; //!< number of diagonals with values

    hmemo::HArray<IndexType> mOffset; //!< size == mNumDiagonals, diagonal offsets
    hmemo::HArray<ValueType> mValues; //!< non-zero values, size is mNumValuesPerDiagonal * mNumDiagonals

private:

    virtual bool checkDiagonalProperty() const;

    // values might be stored row-wise or diagonal-wise

    // take this one for row-major order (fits best for cache architecture)

    // inline IndexType index(IndexType irow, IndexType idiag) const { return irow * mNumDiagonals + idiag; }

    // take this one for diag-major order (fits best for GPUs)

    static inline IndexType diaindex(
        const IndexType irow,
        const IndexType idiag,
        const IndexType numRows,
        const IndexType numDiagonals )
    {
        SCAI_ASSERT_ERROR( irow >= 0, "irow = " << irow );
        SCAI_ASSERT_ERROR( idiag >= 0, "idiag = " << idiag );
        SCAI_ASSERT_ERROR( irow < numRows, "irow = " << irow << " out of range, numRows = " << numRows );
        SCAI_ASSERT_ERROR( idiag < numDiagonals,
                           "idiag = " << idiag << " out of range, numDiagonals = " << numDiagonals );

        return idiag * numRows + irow;
    }

    // Help routine to set offset from used diagonals

    void setOffsets( const IndexType maxNumDiagonals, const bool upperDiagonalUsed[], const bool lowerDiagonalUsed[] );

    static void setUsedDiagonal( bool upperDiagonalUsed[], bool lowerDiagonalUsed[], IndexType i, IndexType j );

    SCAI_LOG_DECL_STATIC_LOGGER( logger ) //!< logger for this matrix format

public:

    // static create method that will be used to register at MatrixStorage factory

    static _MatrixStorage* create();

    // key for factory

    static MatrixCreateKeyType createValue();
};

} /* end namespace lama */

} /* end namespace scai */
