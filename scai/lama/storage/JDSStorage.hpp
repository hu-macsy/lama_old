/**
 * @file JDSStorage.hpp
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
 * @brief Definition of a structure for a (non-distributed) JDS sparse matrix.
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

/** Storage format for a JDS sparse matrix.
 *
 *  The JDS format has the following data:
 *
 *  - number of jagged diagonals (equal to the length of the longest row)
 *  - number of values in each column
 *  - number of values in each row
 *  - original positions of the permuted rows
 *  - values-array
 *
 *  @tparam ValueType is the value type of the matrix values.
 *
 *  Note: default copy constructor and assignment operator can be used.
 */
template<typename ValueType>
class COMMON_DLL_IMPORTEXPORT JDSStorage:
    public MatrixStorage<ValueType>,
    public _MatrixStorage::Register<JDSStorage<ValueType> >    // register at factory
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
    /*  Constructors / Destructor                                           */
    /* ==================================================================== */

    /** Default constructor, creates a storage of size 0 x 0 */

    JDSStorage( hmemo::ContextPtr ctx = hmemo::Context::getContextPtr() );

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
     *   JDSStorage<ValueType> jds( m, n, ctx );
     *   auto jds = zero<JDSStorage<ValueType>>( m, n, ctx );
     *  \endcode
     */
    JDSStorage( const IndexType numRows, const IndexType numColumns, hmemo::ContextPtr ctx = hmemo::Context::getContextPtr() );

    /** Constructor for JDS storage by corresponding arrays. */

    JDSStorage(
        const IndexType numRows,
        const IndexType numColumns,
        hmemo::HArray<IndexType> dlg,
        hmemo::HArray<IndexType> ilg,
        hmemo::HArray<IndexType> perm,
        hmemo::HArray<IndexType> ja,
        hmemo::HArray<ValueType> values,
        hmemo::ContextPtr ctx = hmemo::Context::getContextPtr() );

    /** Destructor of JDS sparse format */

    virtual ~JDSStorage();

    /* ==================================================================== */
    /*  Copy / Move constructor                                             */
    /* ==================================================================== */

    /** Default copy constructor is overridden */

    JDSStorage( const JDSStorage<ValueType>& other );

    /** Move constructor (noexcept allows use in container classes ) */

    JDSStorage( JDSStorage<ValueType>&& other ) noexcept;

    /* ==================================================================== */
    /*  Dynamic / virtual constructors                                      */
    /* ==================================================================== */

    /** Implementation of MatrixStorage::newMatrixStorage for derived class. */

    virtual JDSStorage* newMatrixStorage( const IndexType numRows, const IndexType numColumns ) const;

    JDSStorage* newMatrixStorage() const
    {
        return newMatrixStorage( getNumRows(), getNumColumns() );
    }

    /** Implementation of MatrixStorage::copy for derived class. */

    virtual JDSStorage* copy() const;

    /* ==================================================================== */
    /*   split up member variables                                          */
    /* ==================================================================== */

    void splitUp(
        IndexType& numRows,
        IndexType& numColumns,
        hmemo::HArray<IndexType>& dlg,
        hmemo::HArray<IndexType>& ilg,
        hmemo::HArray<IndexType>& perm,
        hmemo::HArray<IndexType>& ja,
        hmemo::HArray<ValueType>& values );

    /* ==================================================================== */
    /*   assignment operator=                                               */
    /* ==================================================================== */

    /**
     *  Override default assignment operator.
     */
    JDSStorage<ValueType>& operator=( const JDSStorage<ValueType>& other );

    /**
     *  Move assignment operator, reuses allocated data.
     *
     *  The input argument other becomes a zero matrix after successful completion.
     */
    JDSStorage& operator=( JDSStorage<ValueType>&& other );

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
    void assignJDS( const JDSStorage<OtherValueType>& other );

    /* ==================================================================== */
    /*   Other virtual methods                                              */
    /* ==================================================================== */

    /** Implementation of pure method for _MatrixStorage. */

    virtual void clear();

    /** Test the storage data for inconsistencies.
     *
     *  @throw Exception in case of any inconsistency.
     */
    void check( const char* msg ) const;

    /** @brief Allocate arrays ilg and perm
     */
    void allocate( const IndexType numRows, const IndexType numColumns );

    /** Getter routine for the enum value that stands for this format. */

    virtual Format getFormat() const;

    /** Implementation of pure method of class MatrixStorage. */

    virtual void purge();

    /** Implementation of pure method for JDS matrix. */

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
        mepr::StorageWrapper<JDSStorage, SCAI_NUMERIC_TYPES_HOST_LIST>::
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
        mepr::StorageWrapper<JDSStorage, SCAI_NUMERIC_TYPES_HOST_LIST>::
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

    /**
     * Fill up a JDS storage with the given arrays.
     */
    void setJDSData(
        const IndexType numRows,
        const IndexType numColumns,
        const hmemo::HArray<IndexType>& dlg,
        const hmemo::HArray<IndexType>& ilg,
        const hmemo::HArray<IndexType>& perm,
        const hmemo::HArray<IndexType>& ja,
        const hmemo::_HArray& values );

    /** Implementation of MatrixStorage::matrixTimesVector for JDS */

    virtual void matrixTimesVector(
        hmemo::HArray<ValueType>& result,
        const ValueType alpha,
        const hmemo::HArray<ValueType>& x,
        const ValueType beta,
        const hmemo::HArray<ValueType>& y,
        const common::MatrixOp op ) const;

    /** Implementation of MatrixStorage::matrixTimesVectorAsync for JDS */

    virtual tasking::SyncToken* matrixTimesVectorAsync(
        hmemo::HArray<ValueType>& result,
        const ValueType alpha,
        const hmemo::HArray<ValueType>& x,
        const ValueType beta,
        const hmemo::HArray<ValueType>& y,
        const common::MatrixOp op ) const;

    /** Implementation of MatrixStorage::jacobiIterate for JDS */

    virtual void jacobiIterate(
        hmemo::HArray<ValueType>& solution,
        const hmemo::HArray<ValueType>& oldSolution,
        const hmemo::HArray<ValueType>& rhs,
        const ValueType omega ) const;

    /** Implementation of MatrixStorage::jacobiIterateAsync for JDS */

    virtual tasking::SyncToken* jacobiIterateAsync(
        hmemo::HArray<ValueType>& solution,
        const hmemo::HArray<ValueType>& oldSolution,
        const hmemo::HArray<ValueType>& rhs,
        const ValueType omega ) const;

    /** Implementation of MatrixStorage::jacobiIterateHalo for JDS */

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

    /* Print relevant information about matrix storage format. */

    virtual void writeAt( std::ostream& stream ) const;

    const hmemo::HArray<IndexType>& getJA() const;

    const hmemo::HArray<IndexType>& getPerm() const;

    const hmemo::HArray<IndexType>& getDlg() const;

    const hmemo::HArray<IndexType>& getIlg() const;

    const hmemo::HArray<ValueType>& getValues() const;

    IndexType getNumDiagonals() const;

    /** Override _MatrixStorage::getNumValues with a more efficient routine for JDS. */

    virtual IndexType getNumValues() const;

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

    /** Implementation for MatrixStorage::l1Norm */

    virtual RealType<ValueType> l1Norm() const;

    /** Implementation for MatrixStorage::l2Norm */

    virtual RealType<ValueType> l2Norm() const;

    /** Implementation for MatrixStorage::maxNorm */

    virtual RealType<ValueType> maxNorm() const;

    /** Get a value of the matrix.
     *
     * @param[in] i is the row index, 0 <= i < getNumRows()
     * @param[in] j is the colum index, 0 <= j < getNumColumns()
     *
     * Out-of-range check is enabled for DEBUG version.
     */

    ValueType getValue( const IndexType i, const IndexType j ) const;

    /** Implementation of pure method MatrixStorage<ValueType>::setValue for JDS storage */

    void setValue( const IndexType i, const IndexType j, const ValueType val,
                   const common::BinaryOp op = common::BinaryOp::COPY );

    /** Initiate an asynchronous data transfer to a specified location. */

    void prefetch( const hmemo::ContextPtr context ) const;

    /** Will wait for all outstanding asynchronous data transfers. */

    void wait() const;

    /** Swaps this with other.
     * @param[in,out] other the JDSStorage to swap this with
     */
    void swap( JDSStorage<ValueType>& other );

    virtual size_t getMemoryUsageImpl() const;

    using _MatrixStorage::prefetch;
    using _MatrixStorage::getContextPtr;

    using _MatrixStorage::getNumRows;
    using _MatrixStorage::getNumColumns;

    using MatrixStorage<ValueType>::assign;

protected:

    using MatrixStorage<ValueType>::mRowIndexes;
    using MatrixStorage<ValueType>::mCompressThreshold;

private:

    // Note: implicitly: numDiagonals = mDlg.size(), numValues = mJA.size()

    hmemo::HArray<IndexType> mDlg;  //!< number of values in each column, size is mNumDiagonals
    hmemo::HArray<IndexType> mIlg;  //!< number of values in each row, size is getNumRows()
    hmemo::HArray<IndexType> mPerm; //!< position of each row in original matrix, size is getNumRows()
    hmemo::HArray<IndexType> mJA;   //!< column indices, size is mNumValues

    hmemo::HArray<ValueType> mValues; //!< non-zero values (+ optionally zeros in diagonal), size is mNumValues

    /** matrixTimesVector for synchronous and asynchronous execution */

    virtual tasking::SyncToken* gemv(
        hmemo::HArray<ValueType>& result,
        const ValueType alpha,
        const hmemo::HArray<ValueType>& x,
        const ValueType beta,
        const hmemo::HArray<ValueType>& y,
        const common::MatrixOp op,
        bool async ) const;

    /**
     *  @brief Help routine for stable sort of ilg and setting up perm correctly
     */
    void sortRows();

    void print( std::ostream& ) const;

    static std::string initTypeName();

    SCAI_LOG_DECL_STATIC_LOGGER( logger ); //!< logger for this matrix format

};

} /* end namespace lama */

} /* end namespace scai */
