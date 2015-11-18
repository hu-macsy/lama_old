/**
 * @file CSRStorage.hpp
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
 * @brief Definition of a structure for a (non-distributed) CSR sparse matrix.
 * @author Thomas Brandes
 * @date 27.04.2011
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

/** Storage format for a CSR (compressed storage row) sparse matrix.
 *
 *  The CSR format has the following data:
 *
 *  - ia contains the offset for each row to the following arrays
 *  - ja is the array with the corresponding column indexes
 *  - values is the array with all non-zero values
 *
 *  The arrays ja and values are sorted by rows as the array ia indicates.
 *  For each row the matrix elements should be sorted by column indexes
 *  as this implies usually better cache usage. But it is not mandatory.
 *
 *  If the diagonal property flag is set, the first element in each row
 *  is the diagonal element.
 *
 *  \code
 *     values[ja[ia[i]] !! is the matrix diagonal element at (i,i)
 *  \endcode
 *
 *  Usually each row has at least one non-zero element. In case of distributed
 *  matrixes this is not always the case for the non-local or halo part of the
 *  sparse matrix. For such matrixes it might be possible that many rows are
 *  completely empty. In such cases, an additional array rowIndexes will contain the
 *  non-empty rows. This allows faster traversing.
 *
 *  Generally it is not required that the column indexes of one row are sorted.
 *  Nevertheless some operations can benefit from it and so an internal flag
 *  is set if values of one row are sorted by column indexes.
 *
 *  @tparam ValueType is the value type of the matrix values.
 */
template<typename ValueType>
class COMMON_DLL_IMPORTEXPORT CSRStorage: public CRTPMatrixStorage<CSRStorage<ValueType>,ValueType>
{
public:

    typedef ValueType StorageValueType;

    /** get typename of the matrix storage format. */

    static const char* typeName();

    /** Default constructor, stands for a zero matrix of size 0 x 0. */

    CSRStorage();

    /** Constructor for CSR storage by corresponding arrays.
     *
     */

    CSRStorage(
        const IndexType numRows,
        const IndexType numColumns,
        const IndexType numValues,
        const hmemo::LAMAArray<IndexType>& ia,
        const hmemo::LAMAArray<IndexType>& ja,
        const hmemo::ContextArray& values );

    /** Copy constructor can take any matrix storage. */

    explicit CSRStorage( const _MatrixStorage& other )
    {
        assign( other );
    }

    /** Copy constructor can take any matrix storage or context. */

    CSRStorage( const _MatrixStorage& other, const hmemo::ContextPtr context )
    {
        _MatrixStorage::setContextPtr( context );
        assign( other );
    }

    CSRStorage<ValueType>& operator=( const _MatrixStorage& other )
    {
        assign( other );
        return *this;
    }

    /**
     * @brief Implementation of pure method for _MatrixStorage.
     */

    virtual void clear();

    /** Implementation of pure method of class MatrixStorage. */

    virtual void purge();

    /** Implementation of MatrixStorage::copy for derived class. */

    virtual CSRStorage* copy() const;

    /** Implementation of MatrixStorage::create for derived class. */

    virtual CSRStorage* clone() const;

    /** Destructor of CSR sparse storage. */

    virtual ~CSRStorage();

    /** Test the storage data for inconsistencies.
     *
     *  @param[in] msg is additional message string that should be used to identify calling routine
     *  @throws Exception in case of any inconsistency.
     */
    void check( const char* msg ) const;

    /** Getter routine for the enum value that stands for this format. */

    virtual MatrixStorageFormat getFormat() const;

    /** Implementation of pure method.  */

    virtual void setIdentity( const IndexType size );

    //TODO: just temporary

    void setNumValues( const IndexType numValues );

    /**
     * @brief fills CSR sparse matrix by csr sparse data.
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
        const hmemo::LAMAArray<IndexType>& ia,
        const hmemo::LAMAArray<IndexType>& ja,
        const hmemo::LAMAArray<OtherValueType>& values,
        const hmemo::ContextPtr loc ) __attribute__( ( noinline ) );

    /**
     * @brief fills CSR sparse matrix with csr sparse data without a copy operation
     *
     * @param[in] numRows    number of rows
     * @param[in] numColumns number of columns
     * @param[in] numValues  the number of stored elements in the matrix
     * @param[in] ia         row pointer of the input csr sparse matrix (will be lost after)
     * @param[in] ja         column indexes of the input csr sparse matrix (will be lost after)
     * @param[in] values     the data values of the input csr sparse matrix (will be lost after)
     * @param[in] loc        is the context where filling takes place
     */
    template<typename OtherValueType>
    void setCSRDataSwap(
        const IndexType numRows,
        const IndexType numColumns,
        const IndexType numValues,
        hmemo::LAMAArray<IndexType>& ia,
        hmemo::LAMAArray<IndexType>& ja,
        hmemo::LAMAArray<OtherValueType>& values,
        const hmemo::ContextPtr loc );

    /* Print relevant information about matrix storage format. */

    virtual void writeAt( std::ostream& stream ) const;

    /** Getter routine for member variable IA. */

    hmemo::LAMAArray<IndexType>& getIA();

    /** Getter routine for member variable JA. */

    hmemo::LAMAArray<IndexType>& getJA();

    /** Getter routine for member variable values. */

    hmemo::LAMAArray<ValueType>& getValues();

    /** Getter routine for member variable IA (read-only). */

    const hmemo::LAMAArray<IndexType>& getIA() const;

    /** Getter routine for member variable JA (read-only). */

    const hmemo::LAMAArray<IndexType>& getJA() const;

    /** Getter routine for member variable values (read-only). */

    const hmemo::LAMAArray<ValueType>& getValues() const;

    /** Getter routine for the number of stored values. */

    IndexType getNumValues() const;

    /** This method sorts entries in each row by column indexes.
     *
     *  @param[in] diagonalProperty if true first entry in a row is diagonal element if available
     *
     *  This method does not force diagonal property for the storage as it will not fill up
     *  a diagonal element if it is not available.
     *
     *  Note: matrix multiplication with CUSparse requires sorted rows, diagonalProperty = false
     */
    void sortRows( bool diagonalProperty );

    /** Template method for getting row. */

    template<typename OtherType>
    void getRowImpl( hmemo::LAMAArray<OtherType>& row, const IndexType i ) const;

    /** Typed version of getDiagonal
     *
     * @param[out] diagonal is the typed LAMA array for output
     */
    template<typename OtherValueType>
    void getDiagonalImpl( hmemo::LAMAArray<OtherValueType>& diagonal ) const __attribute__( ( noinline ) );

    /** Typed version of setDiagonal
     *
     * @param[in] diagonal is the typed LAMA array for input
     */
    template<typename OtherValueType>
    void setDiagonalImpl( const hmemo::LAMAArray<OtherValueType>& diagonal ) __attribute__( ( noinline ) );

    /** Implementation of pure method. */

    void setDiagonalImpl( const ValueType value );

    /******************************************************************
     *  Scaling of elements in a matrix                                *
     ******************************************************************/

    /** Template version used for virtual routine scale with known value type. */

    template<typename OtherType>
    void scaleImpl( const hmemo::LAMAArray<OtherType>& values ) __attribute__( ( noinline ) );

    /** Implementation of pure method.  */

    void scaleImpl( const ValueType value );

    /** Get a value of the matrix.
     *
     * @param[in] i is the row index, 0 <= i < mNumRows
     * @param[in] j is the colum index, 0 <= j < mNumRows
     *
     * Out-of-range check is enabled for DEBUG version.
     */

    ValueType getValue( const IndexType i, const IndexType j ) const;

    /** Initiate an asynchronous data transfer to a specified location. */

    virtual void prefetch( const hmemo::ContextPtr location ) const;

    /** Will wait for all outstanding asynchronous data transfers. */

    void wait() const;

    /** Allocation of CSR storage with size numRows x numColumns
     *
     * @param[in] numRows    number of rows
     * @param[in] numColumns number of columns
     */
    void allocate( const IndexType numRows, const IndexType numColumns );

    /**
     * @brief compress compresses this by removing non-diagonal elements that
     *        are smaller than eps.
     *
     * compress compresses this by removing non-diagonal elements that are
     * smaller than eps. compress has linear runtime. compress is used by
     * CSRSparseMatrix Multiplication.
     *
     * @param[in] eps   the threshold to set a element of this to zero and
     *                  remove it. (Default: 0.0)
     */
    void compress( const ValueType eps = 0.0 );

    /** Swap this CSR storage data with another CSR storage.
     *
     *  @param[in,out] other is the CSR storage with with this storage is swapped
     *
     *  Note: swap is only possible for two storages of the same format and same type.
     */

    void swap( CSRStorage<ValueType>& other );

    /**
     * @brief Swap the CSR arrays with new arrays.
     *
     * This routine can be used to build a CSR storage with new values. Other member variables
     * (e.g. mDiagonalProperty, rowIndexes, ... ) will be defined correctly.
     */
    void swap( hmemo::LAMAArray<IndexType>& ia, hmemo::LAMAArray<IndexType>& ja, hmemo::LAMAArray<ValueType>& values );

    /** Implementation for pure method is provided. */

    virtual size_t getMemoryUsageImpl() const;

    /**
     * @brief override method of MatrixStorage with a more efficient solution.
     *
     */
    virtual void assignTranspose( const MatrixStorage<ValueType>& other );

    /**
     * @brief override method of MatrixStorage with a more efficient solution.
     *
     */
    virtual void assign( const _MatrixStorage& other );

    /**
     * @brief override method of MatrixStorage with a more efficient solution
     */
    virtual void copyTo( _MatrixStorage& other ) const;

    /** Redistribution of CSR avoids unnecessary conversions. */

    virtual void redistributeCSR( const CSRStorage<ValueType>& other, const class Redistributor& redistributor );

    /** Override splitHalo with version that avoids unnecessary conversions. */

    virtual void splitHalo(
        MatrixStorage<ValueType>& localData,
        MatrixStorage<ValueType>& haloData,
        Halo& halo,
        const Distribution& colDist,
        const Distribution* rowDist ) const;

    /** General routine to build any kind of CSR storage.
     *
     *  @param[out] ia is the CSR offset array
     *  @param[out] ja is the array with the column indexes (optional)
     *  @param[out] values is the array with the non-zero matrix values (optional)
     *  @param[in]  loc is the Context where conversion should be done
     */

    template<typename OtherValueType>
    void buildCSR(
        hmemo::LAMAArray<IndexType>& ia,
        hmemo::LAMAArray<IndexType>* ja,
        hmemo::LAMAArray<OtherValueType>* values,
        const hmemo::ContextPtr loc ) const;

    /**
     *   This routine builds compressed sparse column format data.
     *
     *   Note: this routine does not support type conversion
     */
     
    virtual void buildCSCData(
        hmemo::LAMAArray<IndexType>& colIA,
        hmemo::LAMAArray<IndexType>& colJA,
        hmemo::LAMAArray<ValueType>& cscValues ) const;

    /** Implementation of MatrixStorage::matrixTimesVector for CSR */

    virtual void matrixTimesVector(
        hmemo::LAMAArray<ValueType>& result,
        const ValueType alpha,
        const hmemo::LAMAArray<ValueType>& x,
        const ValueType beta,
        const hmemo::LAMAArray<ValueType>& y ) const;

    /** Implementation of MatrixStorage::vectorTimesMatrix for CSR */
    /** since 1.0.1 */

    virtual void vectorTimesMatrix(
        hmemo::LAMAArray<ValueType>& result,
        const ValueType alpha,
        const hmemo::LAMAArray<ValueType>& x,
        const ValueType beta,
        const hmemo::LAMAArray<ValueType>& y ) const;

    virtual void matrixTimesVectorN(
        hmemo::LAMAArray<ValueType>& result,
        const IndexType n,
        const ValueType alpha,
        const hmemo::LAMAArray<ValueType>& x,
        const ValueType beta,
        const hmemo::LAMAArray<ValueType>& y ) const;

    /** Implementation of MatrixStorage::matrixTimesVectorAsync for CSR */

    virtual tasking::SyncToken* matrixTimesVectorAsync(
        hmemo::LAMAArray<ValueType>& result,
        const ValueType alpha,
        const hmemo::LAMAArray<ValueType>& x,
        const ValueType beta,
        const hmemo::LAMAArray<ValueType>& y ) const;

    /** Implementation of MatrixStorage::vectorTimesMatrixAsync for CSR */
    /** since 1.0.1 */

    virtual tasking::SyncToken* vectorTimesMatrixAsync(
        hmemo::LAMAArray<ValueType>& result,
        const ValueType alpha,
        const hmemo::LAMAArray<ValueType>& x,
        const ValueType beta,
        const hmemo::LAMAArray<ValueType>& y ) const;

    /** Implementation of MatrixStorage::matrixPlusMatrix for CSR */

    virtual void matrixPlusMatrix(
        const ValueType alpha,
        const MatrixStorage<ValueType>& a,
        const ValueType beta,
        const MatrixStorage<ValueType>& b );

    /** Implementation of MatrixStorage::matrixTimesMatrix for CSR */

    virtual void matrixTimesMatrix(
        const ValueType alpha,
        const MatrixStorage<ValueType>& a,
        const MatrixStorage<ValueType>& b,
        const ValueType beta,
        const MatrixStorage<ValueType>& c );

    /** solution = xxx */

    virtual void jacobiIterate(
        hmemo::LAMAArray<ValueType>& solution,
        const hmemo::LAMAArray<ValueType>& oldSolution,
        const hmemo::LAMAArray<ValueType>& rhs,
        const ValueType omega ) const;

    virtual void jacobiIterateHalo(
        hmemo::LAMAArray<ValueType>& localSolution,
        const MatrixStorage<ValueType>& localStorage,
        const hmemo::LAMAArray<ValueType>& haloOldSolution,
        const ValueType omega ) const;

    /** @since 1.1.0 */

    virtual void jacobiIterateHalo(
        hmemo::LAMAArray<ValueType>& localSolution,
        const hmemo::LAMAArray<ValueType>& localDiagonal,
        const hmemo::LAMAArray<ValueType>& haloOldSolution,
        const ValueType omega ) const;

    /** Implementation for MatrixStorage::l1Norm */

    virtual ValueType l1Norm() const;

    /** Implementation for MatrixStorage::l2Norm */

    virtual ValueType l2Norm() const;

    /** Implementation for MatrixStorage::maxNorm */

    virtual ValueType maxNorm() const;

    /** Implementation for MatrixStorage::maxDiffNorm */

    virtual ValueType maxDiffNorm( const MatrixStorage<ValueType>& other ) const;

    /** Implementation of maxDiffNorm for CSR matrices */

    virtual ValueType maxDiffNormImpl( const CSRStorage<ValueType>& other ) const;

    /** print matrix on cout, helpful for debug. */

    void print() const;

    /** Method that returns an array with number of values for each row.
     *
     * @param[out] numRowValues numRowValues[i] is number of stored values in row i
     *
     */

    void buildSparseRowSizes( hmemo::LAMAArray<IndexType>& numRowValues ) const;

    /**
     * This method extracts sparse data (column indexes and data values) for rows
     *
     * @param[out] ja column indexes for the rows
     * @param[out] values corresponding to the columns
     */

    void buildSparseRowData( hmemo::LAMAArray<IndexType>& ja, hmemo::LAMAArray<ValueType>& values ) const;

    using MatrixStorage<ValueType>::prefetch;
    using MatrixStorage<ValueType>::getContextPtr;

protected:

    using MatrixStorage<ValueType>::mNumRows;
    using MatrixStorage<ValueType>::mNumColumns;
    using MatrixStorage<ValueType>::mDiagonalProperty;
    using MatrixStorage<ValueType>::mRowIndexes;
    using MatrixStorage<ValueType>::mCompressThreshold;

    IndexType mNumValues; //!< number of stored elements

    hmemo::LAMAArray<IndexType> mIa; //!< offsets for ja and data, size is numRows+1
    hmemo::LAMAArray<IndexType> mJa; //!< column indexes, size is mIa[ numRows ]
    hmemo::LAMAArray<ValueType> mValues; //!< non-zero values, size is equal to mJa

private:

    bool mSortedRows; //!< if true, all rows are sorted by column indexes
    /**
     * @brief checks if in each row the diagonal element is stored first.
     *
     */
    virtual bool checkDiagonalProperty() const;

    /** Help routine that computes array with row indexes for non-empty rows.
     *  The array is only built if number of non-zero rows is smaller than
     *  a certain percentage ( mThreshold ).
     */

    void buildRowIndexes();

    /** Logger for this class. */

    SCAI_LOG_DECL_STATIC_LOGGER( logger )

    /** Matrix + Matrix for CSR only */

    void    matrixAddMatrixCSR( const ValueType alpha,
                                const CSRStorage<ValueType>& a,
                                const ValueType beta,
                                const CSRStorage<ValueType>& b,
                                const hmemo::ContextPtr loc );

    /** Matrix times Matrix for CSR only */

    void matrixTimesMatrixCSR( const ValueType alpha,
                               const CSRStorage<ValueType>& a,
                               const CSRStorage<ValueType>& b,
                               const hmemo::ContextPtr loc );

    /** result += alpha (*this) * x, where this storage has sparse rows */

    tasking::SyncToken* sparseGEMV( hmemo::LAMAArray<ValueType>& result,
                                    const ValueType alpha,
                                    const hmemo::LAMAArray<ValueType>& x,
                                    bool async ) const;

    /** result = alpha * (*this) * x  */

    tasking::SyncToken* normalGEMV( hmemo::LAMAArray<ValueType>& result,
                                    const ValueType alpha,
                                    const hmemo::LAMAArray<ValueType>& x,
                                    bool async ) const;

    /** result = alpha * (*this) * x + beta * y */

    tasking::SyncToken* normalGEMV( hmemo::LAMAArray<ValueType>& result,
                                    const ValueType alpha,
                                    const hmemo::LAMAArray<ValueType>& x,
                                    const ValueType beta,
                                    const hmemo::LAMAArray<ValueType>& y,
                                    bool async ) const;

    /** matrixTimesVector for synchronous and asynchronous execution */

    virtual tasking::SyncToken* gemv(
        hmemo::LAMAArray<ValueType>& result,
        const ValueType alpha,
        const hmemo::LAMAArray<ValueType>& x,
        const ValueType beta,
        const hmemo::LAMAArray<ValueType>& y,
        bool async ) const;
};

/* --------------------------------------------------------------------------- */

} /* end namespace lama */

} /* end namespace scai */
