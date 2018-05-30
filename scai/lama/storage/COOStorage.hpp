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
#include <scai/lama/storage/MatrixStorage.hpp>
#include <scai/lama/mepr/StorageWrapper.hpp>

// internal scai libraries
#include <scai/sparsekernel/CSRKernelTrait.hpp>
#include <scai/sparsekernel/COOKernelTrait.hpp>
#include <scai/utilskernel/LAMAKernel.hpp>
#include <scai/utilskernel/HArrayUtils.hpp>
#include <scai/utilskernel/UtilKernelTrait.hpp>

#include <memory>

namespace scai
{

namespace lama
{

/** _Matrix storage format for a COO sparse matrix.
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
 * Note: Copy constructor and operator= are supported also for moving data. 
 */
template<typename ValueType>
class COMMON_DLL_IMPORTEXPORT COOStorage:
    public MatrixStorage<ValueType>,
    public _MatrixStorage::Register<COOStorage<ValueType> >    // register at factory
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
    /*  Constructor / Destructor for COOStorage                             */
    /* ==================================================================== */

    /** Default constructor, creates empty storage of size 0 $x$ 0. */

    COOStorage( hmemo::ContextPtr ctx = hmemo::Context::getContextPtr() );

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
     *   COOStorage<ValueType> coo( m, n, ctx );
     *   auto coo = zero<COOStorage<ValueType>>( m, n, ctx );
     *  \endcode
     */
    COOStorage( const IndexType numRows, const IndexType numColumns, hmemo::ContextPtr ctx = hmemo::Context::getContextPtr() );

    /** Constructor for COO storage by passing the corresponding arrays.
     *
     *  @param[in] numRows number of rows
     *  @param[in] numColumns number of columns
     *  @param[in] ia row indexes of the storage entries
     *  @param[in] ja column indexes of the storage entries, ja.size() == ia.size()
     *  @param[in] values matrix values, values.size() == ja.size()
     *  @param[in] ctx specifies the context where the storage is allocated/used
     *
     *  Arrays might also be passed via the move operator in which case the allocated
     *  data is taken over by the constructed COO storage.
     */
    COOStorage(
        const IndexType numRows,
        const IndexType numColumns,
        hmemo::HArray<IndexType> ia,
        hmemo::HArray<IndexType> ja,
        hmemo::HArray<ValueType> values,
        hmemo::ContextPtr ctx = hmemo::Context::getContextPtr() );

    /** Destructor of COO sparse matrix. */

    virtual ~COOStorage();

    /* ==================================================================== */
    /*   copy / move constructors for COOStorage                            */
    /* ==================================================================== */

    /** Default copy constructor is overridden */

    COOStorage( const COOStorage<ValueType>& other );

    /** Default move constructor is overridden */

    COOStorage( COOStorage<ValueType>&& other );

    /* ==================================================================== */
    /*   assignment operator=                                               */
    /* ==================================================================== */

    /**
     *  Override default assignment operator.
     */
    COOStorage<ValueType>& operator=( const COOStorage<ValueType>& other );

    /**
     *  Move assignment operator, reuses allocated data.
     *
     *  The input argument other becomes a zero matrix after successful completion.
     */
    COOStorage& operator=( COOStorage<ValueType>&& other );

    /* ==================================================================== */
    /*   assign methods                                                     */
    /* ==================================================================== */

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
    void assignCOO( const COOStorage<OtherValueType>& other );

    /* ==================================================================== */
    /*   Other methods                                                      */
    /* ==================================================================== */

    virtual void clear();

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
        mepr::StorageWrapper<COOStorage, SCAI_NUMERIC_TYPES_HOST_LIST>::
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
        mepr::StorageWrapper<COOStorage, SCAI_NUMERIC_TYPES_HOST_LIST>::
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

    /* ==================================================================== */
    /*   split up member variables                                         */
    /* ==================================================================== */

    /**
     *  @brief Split this storage up into its member variables
     *
     *  @param numRows    same as getNumRows()
     *  @param numColumns same as getNumColumns()
     *  @param ia         will contain all row indexes
     *  @param ja         will contain all columns indexes
     *  @param values     will contain all non-zero values
     *
     *  Note: this storage is reset to zero by this operation.
     */
    void splitUp(
        IndexType& numRows,
        IndexType& numColumns,
        hmemo::HArray<IndexType>& ia,
        hmemo::HArray<IndexType>& ja,
        hmemo::HArray<ValueType>& values );

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

    /******************************************************************/
    /*  virtual call of default copy constructor                      */
    /******************************************************************/

    /** Implementation of MatrixStorage::copy for derived class. */

    virtual COOStorage* copy() const;

    /** Implementation of MatrixStorage::newMatrixStorage for derived class. */

    virtual COOStorage* newMatrixStorage( const IndexType numRows, const IndexType numColumns ) const;

    virtual COOStorage* newMatrixStorage() const
    {
        return newMatrixStorage( getNumRows(), getNumColumns() );
    }

    /******************************************************************
     *  _Matrix times Vector                                            *
     ******************************************************************/

    /** Implementation of MatrixStorage::matrixTimesVector for COO */

    virtual void matrixTimesVector(
        hmemo::HArray<ValueType>& result,
        const ValueType alpha,
        const hmemo::HArray<ValueType>& x,
        const ValueType beta,
        const hmemo::HArray<ValueType>& y,
        const common::MatrixOp op ) const;

    /** Implementation of MatrixStorage::matrixTimesVectorAsync for COO */

    virtual tasking::SyncToken* matrixTimesVectorAsync(
        hmemo::HArray<ValueType>& result,
        const ValueType alpha,
        const hmemo::HArray<ValueType>& x,
        const ValueType beta,
        const hmemo::HArray<ValueType>& y,
        const common::MatrixOp op ) const;

    /** solution = xxx */

    virtual void jacobiIterate(
        hmemo::HArray<ValueType>& solution,
        const hmemo::HArray<ValueType>& oldSolution,
        const hmemo::HArray<ValueType>& rhs,
        const ValueType omega ) const;

    /** @brief Override default implementation of MatrixStorage::matrixPlusMatrix 
     *
     *  If the target array has COO format, matrix addition is done completely in COO
     *  format. 
     */
    virtual void matrixPlusMatrix(
        const ValueType alpha,
        const MatrixStorage<ValueType>& a,
        const ValueType beta,
        const MatrixStorage<ValueType>& b );

    /** @brief matrixAddMatrix but for COO storages only */

    void matrixPlusMatrixImpl( const ValueType alpha,
                               const COOStorage<ValueType>& a,
                               const ValueType beta,
                               const COOStorage<ValueType>& b );

    /******************************************************************
     *  Halo related operations                                        *
     ******************************************************************/

    /** 
     *  @brief Override MatrixStorage<ValueType>::globalizeHaloIndexes 
     *
     *  This solution is more efficient as temporary CSR data is completely avoided.
     */
    virtual void globalizeHaloIndexes( const dmemo::Halo& halo, const IndexType globalNumColumns );

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
     * @param[in] j is the colum index, 0 <= j < getNumColumns
     *
     * Out-of-range check is enabled for DEBUG version.
     */

    ValueType getValue( const IndexType i, const IndexType j ) const;

    /** Implementation of pure method MatrixStorage<ValueType>::setValue for COO storage */

    void setValue( const IndexType i, const IndexType j, const ValueType val,
                   const common::BinaryOp op = common::BinaryOp::COPY );

    /** Initiate an asynchronous data transfer to a specified location. */

    void prefetch( const hmemo::ContextPtr location ) const;

    /** Will wait for all outstanding asynchronous data transfers. */

    void wait() const;

    /** 
     * @brief Swap this COO storage with another one, i.e. swap all member variables
     *
     * @param[in,out] other the COOStorage to swap this with
     */
    void swap( COOStorage<ValueType>& other );

    virtual size_t getMemoryUsageImpl() const;

    using _MatrixStorage::getNumRows;
    using _MatrixStorage::getNumColumns;
    using _MatrixStorage::getValueType;

    using MatrixStorage<ValueType>::assign;
    using MatrixStorage<ValueType>::prefetch;
    using MatrixStorage<ValueType>::getContextPtr;
    using MatrixStorage<ValueType>::setContextPtr;

protected:

    using MatrixStorage<ValueType>::mRowIndexes;
    using MatrixStorage<ValueType>::mCompressThreshold;

    // coo arrays have all the same size getNumValues()

    hmemo::HArray<IndexType> mIA;     //!< row indices
    hmemo::HArray<IndexType> mJA;     //!< column indices
    hmemo::HArray<ValueType> mValues; //!< non-zero values 

private:

    // check if COO data is sorted; sorts it if not.

    void verifySorting();

    /** matrixTimesVector for synchronous and asynchronous execution */

    virtual tasking::SyncToken* gemv(
        hmemo::HArray<ValueType>& result,
        const ValueType alpha,
        const hmemo::HArray<ValueType>& x,
        const ValueType beta,
        const hmemo::HArray<ValueType>& y,
        const common::MatrixOp op,
        bool async ) const;

    /** result += alpha * (*this) * x */

    tasking::SyncToken* incGEMV( hmemo::HArray<ValueType>& result,
                                 const ValueType alpha,
                                 const hmemo::HArray<ValueType>& x,
                                 const common::MatrixOp op,
                                 bool async ) const;

    SCAI_LOG_DECL_STATIC_LOGGER( logger )

    static std::string initTypeName();

};

} /* end namespace lama */

} /* end namespace scai */
