/**
 * @file COOStorage.hpp
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
#include <scai/utilskernel/HArrayUtils.hpp>
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
    public CRTPMatrixStorage<COOStorage<ValueType>, ValueType>,
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

    /** Template version of getRow */

    template<typename OtherType>
    void getRowImpl( hmemo::HArray<OtherType>& row, const IndexType i ) const;

    /** Template version of setRow */

    template<typename OtherType>
    void setRowImpl( const hmemo::HArray<OtherType>& row, const IndexType i,
                     const utilskernel::binary::BinaryOp op );

    /** Template version of getColumn */

    template<typename OtherType>
    void getColumnImpl( hmemo::HArray<OtherType>& column, const IndexType j ) const;

    /** Template version of setColumn */

    template<typename OtherType>
    void setColumnImpl( const hmemo::HArray<OtherType>& column, const IndexType j,
                        const utilskernel::binary::BinaryOp op );

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

    /** Implementation of pure method MatrixStorage<ValueType>::setValue for ELL storage */

    void setValue( const IndexType i, const IndexType j, const ValueType val,
                   const utilskernel::binary::BinaryOp op = utilskernel::binary::COPY );

    /** Initiate an asynchronous data transfer to a specified location. */

    void prefetch( const hmemo::ContextPtr location ) const;

    /** Will wait for all outstanding asynchronous data transfers. */

    void wait() const;

    /** Swaps this with other.
     * @param[in,out] other the COOStorage to swap this with
     */
    void swapImpl( COOStorage<ValueType>& other );

    /** Implementation for _MatrixStorage::swap */

    virtual void swap( _MatrixStorage& other );

    /**
     * @brief Swap the COO arrays with new arrays.
     *
     * This routine is helpful to get temporay write access to COO arrays. 
     *
     * This routine can be used to build a COO storage with new values. Other member variables
     * will be defined correctly.
     *
     * It is also useful to work on with COO data when COO storage is no more needed.
     */
    void swap( hmemo::HArray<IndexType>& ia, hmemo::HArray<IndexType>& ja, hmemo::HArray<ValueType>& values );

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

} /* end namespace lama */

} /* end namespace scai */
