/**
 * @file CSRStorage.hpp
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
 * @brief Definition of a structure for a (non-distributed) CSR sparse matrix.
 * @author Thomas Brandes
 * @date 27.04.2011
 */

#pragma once

// for dll_import
#include <scai/common/config.hpp>

// base classes
#include <scai/lama/storage/CRTPMatrixStorage.hpp>

// scai libraries
#include <scai/utilskernel/LArray.hpp>
#include <scai/dmemo/Redistributor.hpp>

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
class COMMON_DLL_IMPORTEXPORT CSRStorage:
    public CRTPMatrixStorage<CSRStorage<ValueType>, ValueType>,
    public _MatrixStorage::Register<CSRStorage<ValueType> >    // register at factory
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
        const hmemo::HArray<IndexType>& ia,
        const hmemo::HArray<IndexType>& ja,
        const hmemo::_HArray& values );

    /** Copy constructor that can take any matrix storage. */

    explicit CSRStorage( const _MatrixStorage& other )
    {
        _MatrixStorage::setContextPtr( other.getContextPtr() );
        assign( other );
    }

    /** Override default copy constructor to guarantee same behavior */

    CSRStorage( const CSRStorage<ValueType>& other )

        : CRTPMatrixStorage<CSRStorage<ValueType>, ValueType>()
    {
        _MatrixStorage::setContextPtr( other.getContextPtr() );
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
        // the assignment operator does not change the context any more

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

    /** Implementation of MatrixStorage::newMatrixStorage for derived class. */

    virtual CSRStorage* newMatrixStorage() const;

    /** Destructor of CSR sparse storage. */

    virtual ~CSRStorage();

    /** Test the storage data for inconsistencies.
     *
     *  @param[in] msg is additional message string that should be used to identify calling routine
     *  @throws Exception in case of any inconsistency.
     */
    void check( const char* msg ) const;

    /** Getter routine for the enum value that stands for this format. */

    virtual Format getFormat() const;

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
        const hmemo::HArray<IndexType>& ia,
        const hmemo::HArray<IndexType>& ja,
        const hmemo::HArray<OtherValueType>& values,
        const hmemo::ContextPtr loc ) __attribute__( ( noinline ) );

    /**
     * @brief fills CSR sparse matrix by dia sparse data.
     *
     * @param[in] numRows      number of rows
     * @param[in] numColumns   number of columns
     * @param[in] numDiagonals the number of stored diagonals
     * @param[in] offsets      raw pointer of the input csr sparse matrix
     * @param[in] values       the data values of the input csr sparse matrix
     * @param[in] loc          is the context where filling takes place
     */
    template<typename OtherValueType>
    void setDIADataImpl(
        const IndexType numRows,
        const IndexType numColumns,
        const IndexType numDiagonals,
        const hmemo::HArray<IndexType>& offsets,
        const hmemo::HArray<OtherValueType>& values,
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
        hmemo::HArray<IndexType>& ia,
        hmemo::HArray<IndexType>& ja,
        hmemo::HArray<OtherValueType>& values,
        const hmemo::ContextPtr loc );

    /* Print relevant information about matrix storage format. */

    virtual void writeAt( std::ostream& stream ) const;

    /** Getter routine for member variable IA. */

    utilskernel::LArray<IndexType>& getIA();

    /** Getter routine for member variable JA. */

    utilskernel::LArray<IndexType>& getJA();

    /** Getter routine for member variable values. */

    utilskernel::LArray<ValueType>& getValues();

    /** Getter routine for member variable IA (read-only). */

    const utilskernel::LArray<IndexType>& getIA() const;

    /** Getter routine for member variable JA (read-only). */

    const utilskernel::LArray<IndexType>& getJA() const;

    /** Getter routine for member variable values (read-only). */

    const utilskernel::LArray<ValueType>& getValues() const;

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

    /** This method overrides _MatrixStorage::setDiagonalProperty
     *
     *  This routine only moves the diagonal elements at the beginning of each rows.
     *  It throws an exception if there is no entry for the diagonal element.
     */
    virtual void setDiagonalProperty();

    /** Implementation of pure method MatrixStorage::getSparseRow */

    virtual void getSparseRow( hmemo::HArray<IndexType>& jA, hmemo::_HArray& values, const IndexType i ) const;

    /** Implementation of pure method MatrixStorage::getSparseColumn */

    virtual void getSparseColumn( hmemo::HArray<IndexType>& iA, hmemo::_HArray& values, const IndexType j ) const;

    /** Template version of getRow */

    template<typename OtherType>
    void getRowImpl( hmemo::HArray<OtherType>& row, const IndexType i ) const;

    /** Template version of setRow */

    template<typename OtherType>
    void setRowImpl( const hmemo::HArray<OtherType>& row, const IndexType i,
                     const common::BinaryOp op );

    /** Implementation of pure method MatrixStorage::getColumn */

    void getColumn( hmemo::_HArray& column, const IndexType j ) const;

    /** Template version of setColumn */

    template<typename OtherType>
    void setColumnImpl( const hmemo::HArray<OtherType>& column, const IndexType j,
                        const common::BinaryOp op );

    /** Typed version of getDiagonal
     *
     * @param[out] diagonal is the typed LAMA array for output
     */
    template<typename OtherValueType>
    void getDiagonalImpl( hmemo::HArray<OtherValueType>& diagonal ) const __attribute__( ( noinline ) );

    /** Typed version of setDiagonal
     *
     * @param[in] diagonal is the typed LAMA array for input
     */
    template<typename OtherValueType>
    void setDiagonalImpl( const hmemo::HArray<OtherValueType>& diagonal ) __attribute__( ( noinline ) );

    /** Implementation of pure method. */

    void setDiagonalImpl( const ValueType value );

    /******************************************************************
     *  Scaling of elements in a matrix                                *
     ******************************************************************/

    /** Template version used for virtual routine scale with known value type. */

    template<typename OtherType>
    void scaleImpl( const hmemo::HArray<OtherType>& values ) __attribute__( ( noinline ) );

    /** Implementation of pure method.  */

    void scaleImpl( const ValueType value );

    /** Implementation of pure method.  */

    void conj();

    /** Get a value of the matrix.
     *
     * @param[in] i is the row index, 0 <= i < mNumRows
     * @param[in] j is the colum index, 0 <= j < mNumRows
     *
     * Out-of-range check is enabled for DEBUG version.
     */

    ValueType getValue( const IndexType i, const IndexType j ) const;

    /** Implementation of pure method MatrixStorage<ValueType>::setValue for CSR storage */

    void setValue( const IndexType i, const IndexType j, const ValueType val,
                   const common::BinaryOp op = common::BinaryOp::COPY );

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

    static void compress( 
        hmemo::HArray<IndexType>& ia, 
        hmemo::HArray<IndexType>& ja, 
        hmemo::HArray<ValueType>& values,
        const bool diagonalFlag, 
        const ValueType eps, 
        hmemo::ContextPtr prefContext );

    /** Swap this CSR storage data with another CSR storage.
     *
     *  @param[in,out] other is the CSR storage with with this storage is swapped
     *
     *  Note: swap is only possible for two storages of the same format and same type.
     */

    void swapImpl( CSRStorage<ValueType>& other );

    void swap( _MatrixStorage& other );

    /**
     * @brief Swap the CSR arrays with new arrays.
     *
     * This routine can be used to build a CSR storage with new values. Other member variables
     * (e.g. mDiagonalProperty, rowIndexes, ... ) will be defined correctly.
     */
    void swap( hmemo::HArray<IndexType>& ia, hmemo::HArray<IndexType>& ja, hmemo::HArray<ValueType>& values );

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

    /**
     * @brief override MatrixStorage<ValueType>::copyBlockTo with a more efficient solution
     */
    virtual void copyBlockTo( _MatrixStorage& other, const IndexType first, const IndexType n ) const;

    /** Redistribution of CSR avoids unnecessary conversions. */

    virtual void redistributeCSR( const CSRStorage<ValueType>& other, const dmemo::Redistributor& redistributor );

    /** Override splitHalo with version that avoids unnecessary conversions. */

    virtual void splitHalo(
        MatrixStorage<ValueType>& localData,
        MatrixStorage<ValueType>& haloData,
        dmemo::Halo& halo,
        const dmemo::Distribution& colDist,
        const dmemo::Distribution* rowDist ) const;

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

    /** Own implementation to get global owners of the row */

    virtual void getFirstColumnIndexes( hmemo::HArray<IndexType>& colIndexes ) const;

    /**
     *   This routine builds compressed sparse column format data.
     *
     *   Note: this routine does not support type conversion
     */

    virtual void buildCSCData(
        hmemo::HArray<IndexType>& colIA,
        hmemo::HArray<IndexType>& colJA,
        hmemo::HArray<ValueType>& cscValues ) const;

    /** Implementation of MatrixStorage::matrixTimesVector for CSR */

    virtual void matrixTimesVector(
        hmemo::HArray<ValueType>& result,
        const ValueType alpha,
        const hmemo::HArray<ValueType>& x,
        const ValueType beta,
        const hmemo::HArray<ValueType>& y ) const;

    /** Implementation of MatrixStorage::vectorTimesMatrix for CSR */
    /** since 1.0.1 */

    virtual void vectorTimesMatrix(
        hmemo::HArray<ValueType>& result,
        const ValueType alpha,
        const hmemo::HArray<ValueType>& x,
        const ValueType beta,
        const hmemo::HArray<ValueType>& y ) const;

    virtual void matrixTimesVectorN(
        hmemo::HArray<ValueType>& result,
        const IndexType n,
        const ValueType alpha,
        const hmemo::HArray<ValueType>& x,
        const ValueType beta,
        const hmemo::HArray<ValueType>& y ) const;

    /** Implementation of MatrixStorage::matrixTimesVectorAsync for CSR */

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

    /** binary operation between two CSR storages
     *
     *  @param a, b are the two CSR storages to which the operation is applied to
     *  @param op specifies the binary operator this is applied to entries at same position
     *
     *  Note: op == common::BinaryOp::ADD stands for summing elements up, 
     *        op == common::BinaryOp::COPY for replacing
     *        op == common::BinaryOp::MIN  for minimal values
     */
    void binaryOpCSR( const CSRStorage<ValueType>& a, const CSRStorage<ValueType>& b, common::BinaryOp op );

    /** solution = xxx */

    virtual void jacobiIterate(
        hmemo::HArray<ValueType>& solution,
        const hmemo::HArray<ValueType>& oldSolution,
        const hmemo::HArray<ValueType>& rhs,
        const ValueType omega ) const;

    virtual void jacobiIterateHalo(
        hmemo::HArray<ValueType>& localSolution,
        const MatrixStorage<ValueType>& localStorage,
        const hmemo::HArray<ValueType>& haloOldSolution,
        const ValueType omega ) const;

    /** @since 1.1.0 */

    virtual void jacobiIterateHalo(
        hmemo::HArray<ValueType>& localSolution,
        const hmemo::HArray<ValueType>& localDiagonal,
        const hmemo::HArray<ValueType>& haloOldSolution,
        const ValueType omega ) const;

    /** Implementation for MatrixStorage::l1Norm */

    virtual RealType<ValueType> l1Norm() const;

    /** Implementation for MatrixStorage::l2Norm */

    virtual RealType<ValueType> l2Norm() const;

    /** Implementation for MatrixStorage::maxNorm */

    virtual RealType<ValueType> maxNorm() const;

    /** Implementation for MatrixStorage::maxDiffNorm */

    virtual RealType<ValueType> maxDiffNorm( const MatrixStorage<ValueType>& other ) const;

    /** Implementation of maxDiffNorm for CSR matrices */

    virtual RealType<ValueType> maxDiffNormImpl( const CSRStorage<ValueType>& other ) const;

    /** Implemenation of pure method of class MatrixStorage. */

    virtual void print( std::ostream& ) const;

    /** Method that returns an array with number of values for each row.
     *
     * @param[out] numRowValues numRowValues[i] is number of stored values in row i
     *
     */

    void buildSparseRowSizes( hmemo::HArray<IndexType>& numRowValues ) const;

    /**
     * This method extracts sparse data (column indexes and data values) for rows
     *
     * @param[out] ja column indexes for the rows
     * @param[out] values corresponding to the columns
     */

    void buildSparseRowData( hmemo::HArray<IndexType>& ja, hmemo::HArray<ValueType>& values ) const;

    using MatrixStorage<ValueType>::prefetch;
    using MatrixStorage<ValueType>::getContextPtr;

protected:

    using MatrixStorage<ValueType>::mNumRows;
    using MatrixStorage<ValueType>::mNumColumns;
    using MatrixStorage<ValueType>::mDiagonalProperty;
    using MatrixStorage<ValueType>::mRowIndexes;
    using MatrixStorage<ValueType>::mCompressThreshold;

    IndexType mNumValues; //!< number of stored elements

    utilskernel::LArray<IndexType> mIa; //!< offsets for ja and data, size is numRows+1
    utilskernel::LArray<IndexType> mJa; //!< column indexes, size is mIa[ numRows ]
    utilskernel::LArray<ValueType> mValues; //!< non-zero values, size is equal to mJa

private:

    bool mSortedRows; //!< if true, the column indexes in each row are sorted
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

    /** _Matrix times _Matrix for CSR only */

    void    matrixAddMatrixCSR( const ValueType alpha,
                                const CSRStorage<ValueType>& a,
                                const ValueType beta,
                                const CSRStorage<ValueType>& b,
                                const hmemo::ContextPtr loc );

    /** _Matrix times _Matrix for CSR only */

    void matrixTimesMatrixCSR( const ValueType alpha,
                               const CSRStorage<ValueType>& a,
                               const CSRStorage<ValueType>& b,
                               const hmemo::ContextPtr loc );

    /** result += alpha * (*this) * x, where this storage has sparse rows */

    tasking::SyncToken* sparseGEMV( hmemo::HArray<ValueType>& result,
                                    const ValueType alpha,
                                    const hmemo::HArray<ValueType>& x,
                                    bool async ) const;

    /** result += alpha * x * (*this), where this storage has sparse rows */

    tasking::SyncToken* sparseGEVM( hmemo::HArray<ValueType>& result,
                                    const ValueType alpha,
                                    const hmemo::HArray<ValueType>& x,
                                    bool async ) const;

    /** result = alpha * (*this) * x  */

    tasking::SyncToken* normalGEMV( hmemo::HArray<ValueType>& result,
                                    const ValueType alpha,
                                    const hmemo::HArray<ValueType>& x,
                                    bool async ) const;

    /** result = alpha * x * (*this)  */

    tasking::SyncToken* normalGEVM( hmemo::HArray<ValueType>& result,
                                    const ValueType alpha,
                                    const hmemo::HArray<ValueType>& x,
                                    bool async ) const;

    /** result = alpha * (*this) * x + beta * y */

    tasking::SyncToken* normalGEMV( hmemo::HArray<ValueType>& result,
                                    const ValueType alpha,
                                    const hmemo::HArray<ValueType>& x,
                                    const ValueType beta,
                                    const hmemo::HArray<ValueType>& y,
                                    bool async ) const;

    /** result = alpha * x * (*this) + beta * y */

    tasking::SyncToken* normalGEVM( hmemo::HArray<ValueType>& result,
                                    const ValueType alpha,
                                    const hmemo::HArray<ValueType>& x,
                                    const ValueType beta,
                                    const hmemo::HArray<ValueType>& y,
                                    bool async ) const;

    /** matrixTimesVector for synchronous and asynchronous execution */

    virtual tasking::SyncToken* gemv(
        hmemo::HArray<ValueType>& result,
        const ValueType alpha,
        const hmemo::HArray<ValueType>& x,
        const ValueType beta,
        const hmemo::HArray<ValueType>& y,
        bool async ) const;

    /** vectorTimesMatrix for synchronous and asynchronous execution */

    virtual tasking::SyncToken* gevm(
        hmemo::HArray<ValueType>& result,
        const ValueType alpha,
        const hmemo::HArray<ValueType>& x,
        const ValueType beta,
        const hmemo::HArray<ValueType>& y,
        bool async ) const;

    static std::string initTypeName();

public:

    // static create method that will be used to register at MatrixStorage factory

    static _MatrixStorage* create();

    // key for factory

    static MatrixStorageCreateKeyType createValue();
};

/* --------------------------------------------------------------------------- */

} /* end namespace lama */

} /* end namespace scai */
