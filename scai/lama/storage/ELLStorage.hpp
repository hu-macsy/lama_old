/**
 * @file ELLStorage.hpp
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
 * @brief Definition of a structure for a (non-distributed) ELL sparse matrix.
 * @author Lauretta Schubert
 * @date 25.05.2011
 */

#pragma once

// for dll_import
#include <scai/common/config.hpp>

// base classes
#include <scai/lama/storage/MatrixStorage.hpp>
#include <scai/lama/mepr/StorageWrapper.hpp>

namespace scai
{

namespace lama
{

/** Storage format for an ELLPACK sparse matrix.
 *
 *  The ELLPACK format has the following data:
 *
 *  - values is the array with all non-zero values
 *  - ja is the array with the corresponding column indexes
 *  - ia contains the number of non-zeros in a row
 *
 *  In contrary to the CSR sparse format, the number of values per row
 *  is the same for all rows. So an offset array ia is no more needed.
 *
 *  Therefore, the arrays ja and values might contain more values that
 *  are not needed. These unused values will be set with default values
 *  so that using these entries does not harm. The unused column indexes
 *  should contain a legal column value, the unused entries in values
 *  will be set to zero.
 *
 *  @tparam ValueType is the value type of the matrix values.
 */
template<typename ValueType>
class COMMON_DLL_IMPORTEXPORT ELLStorage:
    public MatrixStorage<ValueType>,
    public _MatrixStorage::Register<ELLStorage<ValueType> >    // register at factory
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

    /* ==================================================================== */
    /*  Constrcutor / Destructor                                            */
    /* ==================================================================== */

    /**
     *  Default constructor, zero matrix.
     */
    ELLStorage( hmemo::ContextPtr ctx = hmemo::Context::getContextPtr() );

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
     *   ELLStorage<ValueType> ell( m, n, ctx );
     *   auto ell = zero<ELLStorage<ValueType>>( m, n, ctx );
     *  \endcode
     */
    ELLStorage( const IndexType numRows, const IndexType numColumns, hmemo::ContextPtr ctx = hmemo::Context::getContextPtr() );

    /** Constructor for ELL storage by corresponding arrays.
     *
     * @param[in] numRows           number of rows
     * @param[in] numColumns        number of columns
     * @param[in] numValuesPerRow   the number of stored elements in the matrix
     * @param[in] ia                row pointer of the input matrix
     * @param[in] ja                column indexes of the input matrix
     * @param[in] values            the data values of the input matrix
     * @param[in] ctx               context for the new storage object
     */
    ELLStorage(
        const IndexType numRows,
        const IndexType numColumns,
        const IndexType numValuesPerRow,
        hmemo::HArray<IndexType> ia,
        hmemo::HArray<IndexType> ja,
        hmemo::HArray<ValueType> values,
        const hmemo::ContextPtr ctx = hmemo::Context::getContextPtr() );

    /** Default copy constructor is overridden */

    ELLStorage( const ELLStorage<ValueType>& other );

    /** Move constructor (noexcept allows use in container classes ) */

    ELLStorage( ELLStorage<ValueType>&& other ) noexcept;

    /** Destructor of ELL sparse matrix. */

    virtual ~ELLStorage();

    /* ==================================================================== */
    /*   assignment operator=                                               */
    /* ==================================================================== */

    /**
     *  Override default assignment operator.
     */
    ELLStorage<ValueType>& operator=( const ELLStorage<ValueType>& other );

    /**
     *  Move assignment operator, reuses allocated data.
     *
     *  The input argument other becomes a zero matrix after successful completion.
     */
    ELLStorage& operator=( ELLStorage<ValueType>&& other );

    /**
     * @brief Implemenation of pure method _MatrixStorage::assign 
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
    void assignELL( const ELLStorage<OtherValueType>& other );

    /* ==================================================================== */
    /*  dynamic constrcutor functions via virtual methods                   */
    /* ==================================================================== */

    /** Implementation of _MatrixStorage::newMatrixStorage using covariant return type */

    virtual ELLStorage* newMatrixStorage( const IndexType numRows, const IndexType numColumns ) const;

    virtual ELLStorage* newMatrixStorage() const
    {
        return newMatrixStorage( getNumRows(), getNumColumns() );
    }

    /** Implementation of _MatrixStorage::copy using covariant return type */

    virtual ELLStorage* copy() const;

    /* ==================================================================== */
    /*   split up member variables                                          */
    /* ==================================================================== */

    /** @brief This method allows the reuse of allocated member variables of this storage.
     *
     *  @param[out] numRows number of rows
     *  @param[out] numColumns number of columns
     *  @param[out] numValuesPerRow number of (maximal) values per row
     *  @param[out] ia  the array with sizes for each row
     *  @param[out] ja  the array with column indexes for the non-zero entries
     *  @param[out] values  the array with values for the non-zero entries
     *
     *  Note: this storage itself is a zero storage afterwards.
     */
    void splitUp(
        IndexType& numRows,
        IndexType& numColumns,
        IndexType& numValuesPerRow,
        hmemo::HArray<IndexType>& ia,
        hmemo::HArray<IndexType>& ja,
        hmemo::HArray<ValueType>& values );

    /**
     *  Implementation of pure method of _MatrixStorage::clear
     */
    virtual void clear();

    /** Test the storage data for inconsistencies.
     *
     *  @param[in] msg is a string message to identify caller site
     *  @throws Exception in case of any inconsistency.
     */
    void check( const char* msg ) const;

    /** Specify sizes of the matrixes for a default matrix. */

    void allocate( const IndexType numRows, const IndexType numColumns );

    /** Getter routine for the enum value that stands for this format. */

    virtual Format getFormat() const;

    /** Implementation of pure method of class MatrixStorage. */

    virtual void purge();

    /** Implementation of pure method. */

    virtual void setIdentity( const IndexType size );

    /* ==================================================================== */
    /*  set / get CSR data                                                  */
    /* ==================================================================== */

    /** Implementation of _MatrixStorage::setCSRData for this class.  */

    void setCSRData(
        const IndexType numRows,
        const IndexType numColumns,
        const hmemo::HArray<IndexType>& ia,
        const hmemo::HArray<IndexType>& ja,
        const hmemo::_HArray& values )
    {
        mepr::StorageWrapper<ELLStorage, SCAI_NUMERIC_TYPES_HOST_LIST>::
            setCSRDataImpl( this, numRows, numColumns, ia, ja, values );
    }

    /**
     * @brief template (non-virtual) version of setCSRData with explicit other value type.
     *
     * @param[in] numRows    number of rows
     * @param[in] numColumns number of columns
     * @param[in] ia         row pointer of the input csr sparse matrix
     * @param[in] ja         column indexes of the input csr sparse matrix
     * @param[in] values     the data values of the input csr sparse matrix
     * @param[in] ctx        is the context where filling takes place
     */
    template<typename OtherValueType>
    void setCSRDataImpl(
        const IndexType numRows,
        const IndexType numColumns,
        const hmemo::HArray<IndexType>& ia,
        const hmemo::HArray<IndexType>& ja,
        const hmemo::HArray<OtherValueType>& values );

    /* ==================================================================== */
    /*  build CSR data                                                      */
    /* ==================================================================== */

    /** Implementation for _MatrixStorage::buildCSRSizes */

    void buildCSRSizes( hmemo::HArray<IndexType>& ia ) const
    {
        hmemo::HArray<IndexType>* ja = NULL;
        hmemo::HArray<ValueType>* values = NULL;
        buildCSR( ia, ja, values, this->getContextPtr() );
    }

    /** Implementation for _MatrixStorage::buildCSRData */

    void buildCSRData( hmemo::HArray<IndexType>& csrIA, hmemo::HArray<IndexType>& csrJA, hmemo::_HArray& csrValues ) const
    {
        mepr::StorageWrapper<ELLStorage, SCAI_NUMERIC_TYPES_HOST_LIST>::
            buildCSRDataImpl( this, csrIA, csrJA, csrValues, getContextPtr() );
    }

    /** 
     *  @brief Template (non-virtual) version of building CSR data
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

    void setELLData(
        const IndexType numRows,
        const IndexType numColumns,
        const IndexType numValuesPerRow,
        const hmemo::HArray<IndexType>& ia,
        const hmemo::HArray<IndexType>& ja,
        const hmemo::_HArray& values );

    /* Print relevant information about matrix storage format. */

    virtual void writeAt( std::ostream& stream ) const;

    /** getter for member variables IA, JA, Data, only const reference */

    const hmemo::HArray<IndexType>& getIA() const;

    const hmemo::HArray<IndexType>& getJA() const;

    const hmemo::HArray<ValueType>& getValues() const;

    /** @brief Getter routine for the number of non-zero values 
     *
     *  This method returns the number of non-zero values, i.e. sum( IA )
     *  but keep in mind that the allocated memory is for getNumValuesPerRow() * getNumRows()
     *  entries and that might be much hginger.
     */
    virtual IndexType getNumValues() const;

    /** Getter routine for the number of num values per row*/

    IndexType getNumValuesPerRow() const;

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
    virtual void getDiagonal( hmemo::HArray<ValueType>& diagonal ) const;

    /** 
     * Implementation of pure method MatrixStorage<ValueType>::setDiagonalV
     */
    virtual void setDiagonalV( const hmemo::HArray<ValueType>& diagonal );

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

    void scale( const ValueType value );

    /** Implementation of pure method.  */

    void conj();

    /** Get a value of the matrix.
     *
     * @param[in] i is the row index, 0 <= i < getNumRows()
     * @param[in] j is the colum index, 0 <= j < getNumColumns()
     *
     * Out-of-range check is enabled for DEBUG version.
     */

    ValueType getValue( const IndexType i, const IndexType j ) const;

    /** Implementation of pure method MatrixStorage<ValueType>::setValue for ELL storage */

    void setValue( const IndexType i, const IndexType j, const ValueType val,
                   const common::BinaryOp op = common::BinaryOp::COPY );

    /** Initiate an asynchronous data transfer to a specified context. */

    void prefetch( const hmemo::ContextPtr context ) const;

    /** Will wait for all outstanding asynchronous data transfers. */

    void wait() const;

    /** Implementation of MatrixStorage::matrixTimesVector for ELL */

    virtual void matrixTimesVector(
        hmemo::HArray<ValueType>& result,
        const ValueType alpha,
        const hmemo::HArray<ValueType>& x,
        const ValueType beta,
        const hmemo::HArray<ValueType>& y,
        const common::MatrixOp op ) const;

    /** Implementation of MatrixStorage::matrixTimesVectorAsync for ELL */

    virtual tasking::SyncToken* matrixTimesVectorAsync(
        hmemo::HArray<ValueType>& result,
        const ValueType alpha,
        const hmemo::HArray<ValueType>& x,
        const ValueType beta,
        const hmemo::HArray<ValueType>& y,
        const common::MatrixOp op ) const;

    /** Implementation of MatrixStorage::jacobiIterate for ELL */

    virtual void jacobiIterate(
        hmemo::HArray<ValueType>& solution,
        const hmemo::HArray<ValueType>& oldSolution,
        const hmemo::HArray<ValueType>& rhs,
        const ValueType omega ) const;

    /** Implementation of MatrixStorage::jacobiIterateAsync for ELL */

    virtual tasking::SyncToken* jacobiIterateAsync(
        hmemo::HArray<ValueType>& solution,
        const hmemo::HArray<ValueType>& oldSolution,
        const hmemo::HArray<ValueType>& rhs,
        const ValueType omega ) const;

    /** Implementation of MatrixStorage::jacobiIterateHalo for ELL  */

    virtual void jacobiIterateHalo(
        hmemo::HArray<ValueType>& localSolution,
        const hmemo::HArray<ValueType>& localDiagonal,
        const hmemo::HArray<ValueType>& haloOldSolution,
        const ValueType omega ) const;

    /** 
     *  @brief Override MatrixStorage<ValueType>::globalizeHaloIndexes 
     *
     *  This solution is more efficient as we need only update of the column data
     */
    virtual void globalizeHaloIndexes( const dmemo::Halo& halo, const IndexType globalNumColumns );

    /** Implementation of MatrixStorage::matrixPlusMatrix for ELL */

    virtual void matrixPlusMatrix(
        const ValueType alpha,
        const MatrixStorage<ValueType>& a,
        const ValueType beta,
        const MatrixStorage<ValueType>& b );

    /** Implementation of MatrixStorage::matrixTimesMatrix for ELL */

    virtual void matrixTimesMatrix(
        const ValueType alpha,
        const MatrixStorage<ValueType>& a,
        const MatrixStorage<ValueType>& b,
        const ValueType beta,
        const MatrixStorage<ValueType>& c );

    /** Implementation for MatrixStorage::l1Norm */

    virtual RealType<ValueType> l1Norm() const;

    /** Implementation for MatrixStorage::l2Norm */

    virtual RealType<ValueType> l2Norm() const;

    /** Implementation for MatrixStorage::maxNorm */

    virtual RealType<ValueType> maxNorm() const;

    /**
     * @brief Override default method MatrixStorage<ValueType>::compress with a more efficient one.
     */
    virtual void compress( const RealType<ValueType> eps = 0 );

    /** Swap this ELL storage data with another ELL storage.
     *
     *  @param[in,out] other is the ELL storage with with this storage is swapped
     *
     *  Note: swap is only possible for two storages of the same format and same type.
     */
    void swap( ELLStorage<ValueType>& other );

    virtual size_t getMemoryUsageImpl() const;

    using _MatrixStorage::getNumRows;
    using _MatrixStorage::getNumColumns;
    using _MatrixStorage::getDiagonalSize;
    using _MatrixStorage::getValueType;

    using MatrixStorage<ValueType>::assign;
    using MatrixStorage<ValueType>::prefetch;
    using MatrixStorage<ValueType>::getContextPtr;
    using MatrixStorage<ValueType>::setContextPtr;

    /** Print out the ELL storage on std::out, use only for debug. */

    virtual void print( std::ostream& ) const;

protected:

    using MatrixStorage<ValueType>::mRowIndexes;
    using MatrixStorage<ValueType>::mCompressThreshold;

private:

    IndexType mNumValuesPerRow; //!< number of values in each row

    hmemo::HArray<IndexType> mIA; //!< size is numRows
    hmemo::HArray<IndexType> mJA; //!< size is numRows x numValuesPerRow
    hmemo::HArray<ValueType> mValues; //!< size is numRows x numValuesPerRow

    /** Addressing function for the arrays ia and ja: column-wise */

    inline IndexType ellindex( const IndexType i, const IndexType jj ) const
    {
        return jj * getNumRows() + i;
    }

    // Alternative addressing row-wise
    // inline IndexType ellindex( const IndexType i, const IndexType jj ) const { return i * mNumValuesPerRow + jj; }

    /** This method fills up unused entries in ja and values with default data to allow more efficient gemv */

    void fillValues();

    /** Help routine that computes array with row indexes for non-empty rows.
     *  The array is only built if number of non-zero rows is smaller than
     *  a certain percentage ( mThreshold ).
     */

    void buildRowIndexes();

    /** Logger for this class. */

    SCAI_LOG_DECL_STATIC_LOGGER( logger )

    /** _Matrix + _Matrix for ELL only */
    void    matrixAddMatrixELL( const ValueType alpha,
                                const ELLStorage<ValueType>& a,
                                const ValueType beta,
                                const ELLStorage<ValueType>& b );

    /** _Matrix times _Matrix for ELL only */
    void matrixTimesMatrixELL( const ValueType alpha,
                               const ELLStorage<ValueType>& a,
                               const ELLStorage<ValueType>& b );

    /** result += alpha (*this) * x, where this storage has sparse rows */

    tasking::SyncToken* sparseGEMV( hmemo::HArray<ValueType>& result,
                                    const ValueType alpha,
                                    const hmemo::HArray<ValueType>& x,
                                    const common::MatrixOp op,
                                    bool async ) const;

    /** result = alpha * (*this) * x  */

    tasking::SyncToken* normalGEMV( hmemo::HArray<ValueType>& result,
                                    const ValueType alpha,
                                    const hmemo::HArray<ValueType>& x,
                                    const common::MatrixOp op,
                                    bool async ) const;

    /** result = alpha * (*this) * x + beta * y */

    tasking::SyncToken* normalGEMV( hmemo::HArray<ValueType>& result,
                                    const ValueType alpha,
                                    const hmemo::HArray<ValueType>& x,
                                    const ValueType beta,
                                    const hmemo::HArray<ValueType>& y,
                                    const common::MatrixOp op,
                                    bool async ) const;

    /** matrixTimesVector for synchronous and asynchronous execution */

    virtual tasking::SyncToken* gemv(
        hmemo::HArray<ValueType>& result,
        const ValueType alpha,
        const hmemo::HArray<ValueType>& x,
        const ValueType beta,
        const hmemo::HArray<ValueType>& y,
        const common::MatrixOp op,
        bool async ) const;

    static std::string initTypeName();

};

} /* end namespace lama */

} /* end namespace scai */
