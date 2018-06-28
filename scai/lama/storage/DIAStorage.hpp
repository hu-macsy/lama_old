/**
 * @file DIAStorage.hpp
 *
 * @license
 * Copyright (c) 2009-2018
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the SCAI framework LAMA.
 *
 * LAMA is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * LAMA is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
 * more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with LAMA. If not, see <http://www.gnu.org/licenses/>.
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
#include <scai/lama/storage/MatrixStorage.hpp>
#include <scai/lama/mepr/StorageWrapper.hpp>

#include <scai/logging.hpp>

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
    public MatrixStorage<ValueType>,
    public _MatrixStorage::Register<DIAStorage<ValueType> >    // register at factory
{
public:

    /* ==================================================================== */
    /*  static getter methods and corresponding pure methods                */
    /* ==================================================================== */

    /** Static method that returns a unique name for this storage class */

    static const char* typeName();

    /** Implementation of pure method _MatrixStorage:getTypeName    */

    virtual const char* getTypeName() const;

    /** Statitc method that return the unique key for matrix storage factory. */

    static MatrixStorageCreateKeyType createValue();

    /** Implementation of pure method _MatrixStorage:getCreateValue    */

    virtual MatrixStorageCreateKeyType getCreateValue() const;

    /** Static method to create a new object of this storage type, used by factory. */

    static _MatrixStorage* create();

    /** Default constructor, creates empty storage of size 0 x 0 */

    DIAStorage( hmemo::ContextPtr ctx = hmemo::Context::getContextPtr() );

    /**
     * @brief Create a zero-storage of a certain size 
     *
     * @param[in] numRows    number of rows
     * @param[in] numColumns number of columns
     * @param[in] ctx        context where storage is located, optional
     *
     * Attention: DEPRECATED.
     *
     * Instead of this constructor you should use the free function zero to create a storage.
     *
     *  \code
     *   DIAStorage<ValueType> dia( m, n, ctx );
     *   auto dia = zero<DIAStorage<ValueType>>( m, n, ctx );
     *  \endcode
     */
    DIAStorage( const IndexType numRows, const IndexType numColumns, hmemo::ContextPtr ctx = hmemo::Context::getContextPtr() );

    /** Constructor for DIA storage by corresponding arrays.
     *
     *  @param[in] numRows    number of rows
     *  @param[in] numColumns number of columns
     *  @param[in] offsets    indexes of diagonals, size is number of diagonals
     *  @param[in] values     values of the diagonals, size is offsets.size() * numrows
     *  @param[in] ctx        context for the storage
     *
     *  \code
     *      DIAStorage<double> dia( 3, 3, HArray<int>( { 0, 1 } ), HArray<double>( { 1.0, 0.5, 1.0, 0.5, 1.0, 0.0 } ) );
     *  \endcode
     */
    DIAStorage(
        const IndexType numRows,
        const IndexType numColumns,
        hmemo::HArray<IndexType> offsets,
        hmemo::HArray<ValueType> values,
        hmemo::ContextPtr ctx = hmemo::Context::getContextPtr() );

    /** Default copy constructor is overridden */

    DIAStorage( const DIAStorage<ValueType>& other );

    /** Move constructor (noexcept allows use in container classes ) */

    DIAStorage( DIAStorage<ValueType>&& other ) noexcept;

    /** Implementation of MatrixStorage::copy for derived class. */

    virtual DIAStorage* copy() const;

    /** Implementation of MatrixStorage::newMatrixStorage for derived class. */

    virtual DIAStorage* newMatrixStorage( const IndexType numRows, const IndexType numColumns ) const;

    virtual DIAStorage* newMatrixStorage() const
    {
        return newMatrixStorage( getNumRows(), getNumColumns() );
    }

    /** Implementation of _MatrixStorage::clear  */

    virtual void clear();

    /** Destructor of DIA sparse matrix. */

    virtual ~DIAStorage();

    /* ==================================================================== */
    /*   assignment operator=                                               */
    /* ==================================================================== */

    /**
     *  Override default assignment operator.
     */
    DIAStorage<ValueType>& operator=( const DIAStorage<ValueType>& other );

    /**
     *  Move assignment operator, reuses allocated data.
     *
     *  The input argument other becomes a zero matrix after successful completion.
     */
    DIAStorage& operator=( DIAStorage<ValueType>&& other );

    /**
     * @brief Implementation of pure method _MatrixStorage::assign
     */
    virtual void assign( const _MatrixStorage& other );

    /**
     * @brief Implemenation of pure method MatrixStorage<ValueType>::assignDiagonal
     */
    virtual void assignDiagonal( const hmemo::HArray<ValueType>& diagonal );

    /**
     *  @brief Implemenation of assignments for this class
     */
    template<typename OtherValueType>
    void assignImpl( const MatrixStorage<OtherValueType>& other );

    /**
     *  @brief Implementation of assign method for same storage type. 
     */
    template<typename OtherValueType>
    void assignDIA( const DIAStorage<OtherValueType>& other );


    /* ==================================================================== */
    /*   Implementation of other pure methods                               */
    /* ==================================================================== */

    /** Test the storage data for inconsistencies.
     *
     *  @throw Exception in case of any inconsistency.
     */
    void check( const char* msg ) const;

    /** Getter routine for the enum value that stands for this format. */

    virtual Format getFormat() const;

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

    /* ==================================================================== */
    /*  set / get CSR data                                                  */
    /* ==================================================================== */

    /** Implementation of _MatrixStorage::setCSRData for this class.  
     *
     *  Context for conversion is given by the context of the storage.
     */
    void setCSRData(
        const IndexType numRows,
        const IndexType numColumns,
        const hmemo::HArray<IndexType>& ia,
        const hmemo::HArray<IndexType>& ja,
        const hmemo::_HArray& values );

    /**
     * @brief template (non-virtual) version of setCSRData with explicit other value type.
     *
     * @param[in] numRows    number of rows
     * @param[in] numColumns number of columns
     * @param[in] ia         row pointer of the input csr sparse matrix
     * @param[in] ja         column indexes of the input csr sparse matrix
     * @param[in] values     the data values of the input csr sparse matrix
     */
    void setCSRDataImpl(
        const IndexType numRows,
        const IndexType numColumns,
        const hmemo::HArray<IndexType>& ia,
        const hmemo::HArray<IndexType>& ja,
        const hmemo::HArray<ValueType>& values );

    /* ==================================================================== */
    /*  build CSR data                                                      */
    /* ==================================================================== */

    /** Implementation for _MatrixStorage::buildCSRSizes */

    void buildCSRSizes( hmemo::HArray<IndexType>& ia ) const;

    /** Implementation for _MatrixStorage::buildCSRData */

    void buildCSRData( hmemo::HArray<IndexType>& csrIA, hmemo::HArray<IndexType>& csrJA, hmemo::_HArray& csrValues ) const;

    /** Implementation of MatrixStorage::matrixTimesVector for DIA */

    virtual void matrixTimesVector(
        hmemo::HArray<ValueType>& result,
        const ValueType alpha,
        const hmemo::HArray<ValueType>& x,
        const ValueType beta,
        const hmemo::HArray<ValueType>& y,
        const common::MatrixOp op ) const;

    /** Implementation of MatrixStorage::matrixTimesVectorAsync for DIA */

    virtual tasking::SyncToken* matrixTimesVectorAsync(
        hmemo::HArray<ValueType>& result,
        const ValueType alpha,
        const hmemo::HArray<ValueType>& x,
        const ValueType beta,
        const hmemo::HArray<ValueType>& y,
        const common::MatrixOp op ) const;

    /** Implementation of MatrixStorage::jacobiIterate for DIA */

    virtual void jacobiIterate(
        hmemo::HArray<ValueType>& solution,
        const hmemo::HArray<ValueType>& oldSolution,
        const hmemo::HArray<ValueType>& rhs,
        const ValueType omega ) const;

    /** Implementation of MatrixStorage::jacobiIterateHalo for DIA */

    virtual void jacobiIterateHalo(
        hmemo::HArray<ValueType>& localSolution,
        const hmemo::HArray<ValueType>& localDiagonal,
        const hmemo::HArray<ValueType>& haloOldSolution,
        const ValueType omega ) const;

    /* Print relevant information about matrix storage format. */

    virtual void writeAt( std::ostream& stream ) const;

    /** Getter routine for the matrix offsets (read-only). */

    const hmemo::HArray<IndexType>& getOffsets() const;

    /** Getter routine for the matrix values (read-only). */

    const hmemo::HArray<ValueType>& getValues() const;

    /** Getter routine for the number of diagonals. */

    IndexType getNumDiagonals() const;

    /******************************************************************/
    /*  set - get  row - column                                       */
    /******************************************************************/

    /** Implementation of pure method MatrixStorage<ValueType>::getRow */

    virtual void getRow( hmemo::HArray<ValueType>& row, const IndexType i ) const;

    /** Implementation of pure method MatrixStorage<ValueType>::getColumn */

    virtual void getColumn( hmemo::HArray<ValueType>& column, const IndexType j ) const;

    /** Implementation of pure method MatrixStorage<ValueType>::getSparseRow */

    virtual void getSparseRow( hmemo::HArray<IndexType>& jA, hmemo::HArray<ValueType>& values, const IndexType i ) const;

    /** Implementation of pure method MatrixStorage::getSparseColumn */

    virtual void getSparseColumn( hmemo::HArray<IndexType>& iA, hmemo::HArray<ValueType>& values, const IndexType j ) const;

    /** Implementation of pure method MatrixStorage<ValueType>::setRow */

    virtual void setRow( const hmemo::HArray<ValueType>& row, const IndexType i, const common::BinaryOp op );

    /** Implementation of pure method MatrixStorage<ValueType>::setColumn */

    virtual void setColumn( const hmemo::HArray<ValueType>& column, const IndexType j, const common::BinaryOp op );

    /******************************************************************/
    /*  set / get diagonal                                            */
    /******************************************************************/

    /** 
     * Implementation of pure method MatrixStorage<ValueType>::getDiagonal
     */
    void getDiagonal( hmemo::HArray<ValueType>& diagonal ) const;

    /** 
     * Implementation of pure method MatrixStorage<ValueType>::setDiagonalV
     */
    void setDiagonalV( const hmemo::HArray<ValueType>& diagonal );

    /** 
     * Implementation of pure method MatrixStorage<ValueType>::setDiagonal
     */
    virtual void setDiagonal( const ValueType value );

    /******************************************************************
     *  Scaling of elements in a matrix                                *
     ******************************************************************/

    /** Template version used for virtual routine scale with known value type. */

    void scaleRows( const hmemo::HArray<ValueType>& values );

    /** Implementation of pure method.  */

    virtual void scale( const ValueType value );

    /** Implementation of pure method.  */

    virtual void conj();

    /** Implementation for MatrixStorage::l1Norm */

    virtual RealType<ValueType> l1Norm() const;

    /** Implementation for MatrixStorage::l2Norm */

    virtual RealType<ValueType> l2Norm() const;

    /** Implementation for MatrixStorage::maxNorm */

    virtual RealType<ValueType> maxNorm() const;

    /** Implementation of pure method. */

    ValueType getValue( const IndexType i, const IndexType j ) const;

    /** Implementation of pure method MatrixStorage<ValueType>::setValue for DIA storage */

    void setValue( const IndexType i, const IndexType j, const ValueType val,
                   const common::BinaryOp op = common::BinaryOp::COPY );

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

    void print( std::ostream& stream = std::cout ) const;

    using _MatrixStorage::getNumRows;
    using _MatrixStorage::getNumColumns;
    using _MatrixStorage::getValueType;

    using MatrixStorage<ValueType>::prefetch;
    using MatrixStorage<ValueType>::getContextPtr;
    using MatrixStorage<ValueType>::assign;
    using MatrixStorage<ValueType>::setContextPtr;

protected:

    using MatrixStorage<ValueType>::mRowIndexes;
    using MatrixStorage<ValueType>::mCompressThreshold;

private:

    hmemo::HArray<IndexType> mOffset; //!< offset for each stored diagonal, entries are sorted
    hmemo::HArray<ValueType> mValues; //!< non-zero values, size is getNumRows() * mOffset.size()

    // values might be stored row-wise or diagonal-wise

    // take this one for row-major order (fits best for cache architecture)

    // inline IndexType index(IndexType irow, IndexType idiag) const { return irow * numDiagonals + idiag; }

    // take this one for diag-major order (fits best for GPUs)

    static inline IndexType diaindex(
        const IndexType irow,
        const IndexType idiag,
        const IndexType numRows,
        const IndexType numDiagonals )
    {
        SCAI_ASSERT_VALID_INDEX_ERROR( irow, numRows, "illegal row index" )
        SCAI_ASSERT_VALID_INDEX_ERROR( idiag, numDiagonals, "illegal diag" )

        return idiag * numRows + irow;
    }

    /** matrixTimesVector for synchronous and asynchronous execution */

    virtual tasking::SyncToken* gemv(
        hmemo::HArray<ValueType>& result,
        const ValueType alpha,
        const hmemo::HArray<ValueType>& x,
        const ValueType beta,
        const hmemo::HArray<ValueType>& y,
        const common::MatrixOp op,
        bool async ) const;

    // Help routine to find the main diagonal

    IndexType getMainIndex() const;

    static std::string initTypeName();

    SCAI_LOG_DECL_STATIC_LOGGER( logger ); //!< logger for this matrix format

};

} /* end namespace lama */

} /* end namespace scai */
