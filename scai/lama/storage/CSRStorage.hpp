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
#include <scai/lama/storage/MatrixStorage.hpp>
#include <scai/lama/mepr/StorageWrapper.hpp>

// scai libraries
#include <scai/hmemo/HArray.hpp>
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
    public MatrixStorage<ValueType>,
    public _MatrixStorage::Register<CSRStorage<ValueType> >    // register at factory
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
    /*  Constructor / Desctructor                                           */
    /* ==================================================================== */

    /** Default constructor, stands for a zero matrix of size 0 x 0. 
     *  
     *  @param[in] ctx specifies the context of the storage, optional, defaults to current context
     */
    CSRStorage( hmemo::ContextPtr ctx = hmemo::Context::getContextPtr() );

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
     *   CSRStorage<ValueType> csr( m, n, ctx );
     *   auto csr = zero<CSRStorage<ValueType>>( m, n, ctx );
     *  \endcode
     */
    CSRStorage( const IndexType numRows, const IndexType numColumns, hmemo::ContextPtr ctx = hmemo::Context::getContextPtr() );

    /** Constructor for CSR storage by passing the corresponding arrays.
     *
     *  @param[in] numRows number of rows
     *  @param[in] numColumns number of columns
     *  @param[in] ia is the offset array, size must be numRows + 1
     *  @param[in] ja column indexes
     *  @param[in] values matrix values, values.size() == ja.size()
     *  @param[in] ctx specifies optionally the context
     *
     *  Arrays might also be passed via the move operator in which case the allocated
     *  data is taken over by the constructed CSR storage.
     */
    CSRStorage(
        const IndexType numRows,
        const IndexType numColumns,
        hmemo::HArray<IndexType> ia,
        hmemo::HArray<IndexType> ja,
        hmemo::HArray<ValueType> values,
        hmemo::ContextPtr ctx = hmemo::Context::getContextPtr() );

    /** Destructor of CSR sparse storage. */

    virtual ~CSRStorage();

    /* ==================================================================== */
    /*   Copy / move constructors                                           */
    /* ==================================================================== */

    /** Override default copy constructor to guarantee same behavior */

    CSRStorage( const CSRStorage<ValueType>& other );

    /** Move constructor
     *
     *  The move constructor must be defined with attribute noexcept to allow its use in
     *  the C++ container classes. Exceptions thrown inside (e.g. if there are any 
     *  read and write accesses to the involved arrays) cause program termination.
     */
    CSRStorage( CSRStorage<ValueType>&& other ) noexcept;

    /* ==================================================================== */
    /*   assignment operator=                                               */
    /* ==================================================================== */

    /**
     *  Override default assignment operator.
     */
    CSRStorage<ValueType>& operator=( const CSRStorage<ValueType>& other );

    /**
     *  Move assignment operator, reuses allocated data.
     *
     *  The input argument other becomes a zero matrix after successful completion.
     */
    CSRStorage& operator=( CSRStorage<ValueType>&& other );

    /* ==================================================================== */
    /*  set CSR data                                                        */
    /* ==================================================================== */

    /** Implementation of _MatrixStorage::setCSRData for this class.  */

    void setCSRData(
        const IndexType numRows,
        const IndexType numColumns,
        const hmemo::HArray<IndexType>& ia,
        const hmemo::HArray<IndexType>& ja,
        const hmemo::_HArray& values )
    {
        mepr::StorageWrapper<CSRStorage, SCAI_NUMERIC_TYPES_HOST_LIST>::
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

    void buildCSRSizes( hmemo::HArray<IndexType>& ia ) const;

    /** Implementation for _MatrixStorage::buildCSRData */

    void buildCSRData( hmemo::HArray<IndexType>& csrIA, hmemo::HArray<IndexType>& csrJA, hmemo::_HArray& csrValues ) const;

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
    /*   split up member variables                                          */
    /* ==================================================================== */

    void splitUp(
        IndexType& numRows,
        IndexType& numColumns,
        hmemo::HArray<IndexType>& ia,
        hmemo::HArray<IndexType>& ja,
        hmemo::HArray<ValueType>& values );

    /* ==================================================================== */
    /*  Implementation of pure methods                                      */
    /* ==================================================================== */

    /**
     * @brief Implementation of pure method _MatrixStorage.
     */

    virtual void clear();

    /** Implementation of pure method of class MatrixStorage. */

    virtual void purge();

    /** Implementation of MatrixStorage::copy for derived class. */

    virtual CSRStorage* copy() const;

    /** Implementation of MatrixStorage::newMatrixStorage for derived class. */

    virtual CSRStorage* newMatrixStorage( const IndexType numRows, const IndexType numColumns ) const;

    virtual CSRStorage* newMatrixStorage() const
    {
        return newMatrixStorage( getNumRows(), getNumColumns() );
    }

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

    /* Print relevant information about matrix storage format. */

    virtual void writeAt( std::ostream& stream ) const;

    /** Getter routine for member variable IA (read-only). */

    const hmemo::HArray<IndexType>& getIA() const;

    /** Getter routine for member variable JA (read-only). */

    const hmemo::HArray<IndexType>& getJA() const;

    /** Getter routine for member variable values (read-only). */

    const hmemo::HArray<ValueType>& getValues() const;

    /** Getter routine for the number of stored values. */

    IndexType getNumValues() const;

    /** This method sorts entries in each row by column indexes.
     *
     *  Note: matrix multiplication with CUSparse requires sorted rows
     */
    void sortRows();

    /** 
     *  @brief Query if the entries in each row are sorted according to the column indexes. 
     */
    bool hasSortedRows();

    /** 
     *  @brief Set the diagonal element as first entry in each row.
     * 
     *  This operation might be used to optimize the access to the diagonal of the storage.
     *  The methods setDiagonal(V) and getDiagonal might be significantly faster. 
     */
    void setDiagonalFirst();

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
     *
     * This method might be faster if either the diagonal element is the first entry in a row
     * or if the entries of each row are sorted.
     */
    virtual void getDiagonal( hmemo::HArray<ValueType>& diagonal ) const;

    /** 
     * Implementation of pure method MatrixStorage<ValueType>::setDiagonalV
     *
     * This method does not change the sparsity pattern. So it throws an
     * exception if a non-zero value is set for a non existing entry.
     */
    virtual void setDiagonalV( const hmemo::HArray<ValueType>& diagonal );

    /** 
     * Implementation of pure method MatrixStorage<ValueType>::setDiagonal
     *
     * This method does not change the sparsity pattern. So it throws an
     * exception if a non-zero value is set for a non existing entry.
     */
    virtual void setDiagonal( const ValueType value );

    /*******************************************************************
     *  Scaling of elements in a matrix                                *
     ******************************************************************/

    /** Implementation of pure method MatrixStorage<ValueType>::scaleRows */

    void scaleRows( const hmemo::HArray<ValueType>& values );

    /** Implementation of pure method MatrixStorage<ValueType>::scaleColumns */

    void scaleColumns( const hmemo::HArray<ValueType>& values );

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
     *  @brief Override default implementation MatrixStorage<ValueType>::compress 
     */
    virtual void compress( const RealType<ValueType> eps = 0 );

    /** Swap this CSR storage data with another CSR storage.
     *
     *  @param[in,out] other is the CSR storage with with this storage is swapped
     */
    void swap( CSRStorage<ValueType>& other );

    /**
     * @brief Swap the CSR arrays with new arrays.
     *
     * This routine can be used to build a CSR storage with new values. Other member variables
     * will be defined correctly.
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
    void assignCSR( const CSRStorage<OtherValueType>& other );

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

    /** 
     *  @brief Override MatrixStorage<ValueType>::globalizeHaloIndexes 
     *
     *  This solution is more efficient as temporary CSR data is completely avoided.
     */ 
    virtual void globalizeHaloIndexes( const dmemo::Halo& halo, const IndexType globalNumColumns );

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
        const hmemo::HArray<ValueType>& y,
        const common::MatrixOp op ) const;

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
        const hmemo::HArray<ValueType>& y, 
        const common::MatrixOp op ) const;

    /** Implementation of element-wise binary operation for CSR */

    virtual void binaryOp(
        const MatrixStorage<ValueType>& a,
        const common::BinaryOp op,
        const MatrixStorage<ValueType>& b );

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
    void binaryOpCSR( const CSRStorage<ValueType>& a, common::BinaryOp op, const CSRStorage<ValueType>& b );

    /** solution = xxx */

    virtual void jacobiIterate(
        hmemo::HArray<ValueType>& solution,
        const hmemo::HArray<ValueType>& oldSolution,
        const hmemo::HArray<ValueType>& rhs,
        const ValueType omega ) const;

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

    /* ========================================================================= */
    /*       Filling matrix storage with assembled COO data                      */
    /* ========================================================================= */

    /**
     *  @brief Override MatrixStorage<ValueType>::fillCOO 
     *
     *  This method 'must' override the method of the base class
     */
    virtual void fillCOO(
        hmemo::HArray<IndexType> ia,
        hmemo::HArray<IndexType> ja,
        hmemo::HArray<ValueType> values,
        const common::BinaryOp op = common::BinaryOp::COPY );

    using _MatrixStorage::prefetch;
    using _MatrixStorage::getContextPtr;

    using _MatrixStorage::getNumRows;
    using _MatrixStorage::getNumColumns;

    using MatrixStorage<ValueType>::getValueType;

private:

    /* ==================================================================== */
    /*  Member variables ( always private )                                 */
    /* ==================================================================== */

    hmemo::HArray<IndexType> mIA; //!< offsets for ja and data, size is numRows+1
    hmemo::HArray<IndexType> mJA; //!< column indexes, size is mIA[ numRows ]
    hmemo::HArray<ValueType> mValues; //!< non-zero values, size is equal to mJA

    bool mSortedRows; //!< if true, the column indexes in each row are sorted

    using MatrixStorage<ValueType>::mRowIndexes;
    using MatrixStorage<ValueType>::mCompressThreshold;

    static std::string initTypeName();

    /** Help routine that computes array with row indexes for non-empty rows.
     *  The array is only built if number of non-zero rows is smaller than
     *  a certain percentage ( mThreshold ).
     */
    void buildRowIndexes();

    /** Logger for this class. */

    SCAI_LOG_DECL_STATIC_LOGGER( logger )

    /** @brief matrixPlusMatrix but for CSR storages only, no more virtual */

    void matrixPlusMatrixImpl( const ValueType alpha,
                               const CSRStorage<ValueType>& a,
                               const ValueType beta,
                               const CSRStorage<ValueType>& b );

    /** _Matrix times _Matrix for CSR only */

    void matrixTimesMatrixCSR( const ValueType alpha,
                               const CSRStorage<ValueType>& a,
                               const CSRStorage<ValueType>& b );

    /** matrixTimesVector for synchronous and asynchronous execution */

    virtual tasking::SyncToken* gemv(
        hmemo::HArray<ValueType>& result,
        const ValueType alpha,
        const hmemo::HArray<ValueType>& x,
        const ValueType beta,
        const hmemo::HArray<ValueType>& y,
        const common::MatrixOp op,
        bool async ) const;
};

/* --------------------------------------------------------------------------- */
/*   Implementation of inline methods                                          */
/* --------------------------------------------------------------------------- */

template<typename ValueType>
CSRStorage<ValueType>::CSRStorage( const CSRStorage<ValueType>& other ) :

    MatrixStorage<ValueType>( other )

{
    mIA     = other.mIA;
    mJA     = other.mJA;
    mValues = other.mValues;

    mSortedRows       = other.mSortedRows;
    buildRowIndexes();
}

template<typename ValueType>
CSRStorage<ValueType>::CSRStorage( CSRStorage<ValueType>&& other ) noexcept :

    MatrixStorage<ValueType>( std::move( other ) ),
    mIA( std::move( other.mIA ) ),
    mJA( std::move( other.mJA ) ),
    mValues( std::move( other.mValues ) )
{
    mSortedRows = other.mSortedRows;
}

/* --------------------------------------------------------------------------- */

} /* end namespace lama */

} /* end namespace scai */
