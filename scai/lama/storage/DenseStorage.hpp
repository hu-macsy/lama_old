/**
 * @file DenseStorage.hpp
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
 * @brief Definition of a structure for a dense matrix.
 * @author Thomas Brandes, Lauretta Schubert
 * @date 27.04.2011
 */

#pragma once

// for dll_import
#include <scai/common/config.hpp>

// base classes
#include <scai/lama/storage/MatrixStorage.hpp>
#include <scai/lama/mepr/StorageWrapper.hpp>

// internal scai libraries
#include <scai/hmemo.hpp>

namespace scai
{

namespace lama
{

template<typename ValueType> class DenseStorage;

/**
 * Storage format for a Dense matrix.
 *
 *  @tparam ValueType is the value type of the matrix values.
 */
template<typename ValueType>
class COMMON_DLL_IMPORTEXPORT DenseStorage:
    public MatrixStorage<ValueType>,
    public _MatrixStorage::Register<DenseStorage<ValueType> >    // register at factory
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
    /*   Constructor / Destructor                                           */
    /* ==================================================================== */

    /** @brief Default constructor for a dense storage of 0 x 0  */

    DenseStorage( hmemo::ContextPtr ctx = hmemo::Context::getContextPtr() );

    /**
     * @brief Create a storage of a certain size 
     *
     * @param[in] numRows    number of rows
     * @param[in] numColumns number of columns
     * @param[in] ctx        specifies optionally the context where storage is located
     *
     * Attention: DEPRECATED
     *
     *  \code
     *   DenseStorage<ValueType> dense( m, n, ctx );
     *   auto dense = undefined<DenseStorage<ValueType>>( m, n, ctx );
     *   DenseStorage<ValueType> dense( m, n, ctx );
     *   dense.setZero();
     *   auto dense = zero<DenseStorage<ValueType>>( m, n, ctx );
     *  \endcode
     */
    DenseStorage( const IndexType numRows, const IndexType numColumns, hmemo::ContextPtr ctx = hmemo::Context::getContextPtr() );

    /** Constructor for dense storage by passing the shape and the 'dense' array with values 
     *
     *  @param[in] numRows number of rows
     *  @param[in] numColumns number of columns
     *  @param[in] data contains linearized matrix values (row major order), data.size() = numRows * numColumns
     *  @param[in] ctx  is the context of the new constructed storage
     *
     *  The data array might also be passed via the move operator in which case the allocated
     *  data is taken over by the constructed dense storage.
     */
    DenseStorage( 
        const IndexType numRows, 
        const IndexType numColumns, 
        hmemo::HArray<ValueType> data,
        hmemo::ContextPtr ctx = hmemo::Context::getContextPtr() );

    /** Default copy constructor is overridden */

    DenseStorage( const DenseStorage<ValueType>& other );

    /** Default move constructor is overridden */

    DenseStorage( DenseStorage<ValueType>&& other );

    /** Destructor of Dense matrix reference. */

    virtual ~DenseStorage();

    /* ==================================================================== */
    /*   assignment operator=                                               */
    /* ==================================================================== */

    /**
     *  Override default assignment operator.
     */
    DenseStorage<ValueType>& operator=( const DenseStorage<ValueType>& other );

    /**
     *  Move assignment operator, reuses allocated data.
     *
     *  The input argument other becomes a zero matrix after successful completion.
     */
    DenseStorage& operator=( DenseStorage<ValueType>&& other );

    /** Implementation of pure method _MatrixStorage::assign */

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
    void assignDense( const DenseStorage<OtherValueType>& other );

    /* ==================================================================== */
    /*   Split up for reuse of member variables                             */
    /* ==================================================================== */

    /**
     *  @brief extract data values from dense storage 
     *
     *  This dense storage itself behaves afterwards in the same way as it has been moved.
     */
    void splitUp( IndexType& numRows, IndexType& numColumns, hmemo::HArray<ValueType>& data );

    /* ==================================================================== */
    /*   Dynamic (virtual) constructors                                     */
    /* ==================================================================== */

    /** Implementation of MatrixStorage::newMatrixStorage for derived class. */

    virtual DenseStorage<ValueType>* newMatrixStorage( const IndexType numRows, const IndexType numColumns ) const;

    DenseStorage* newMatrixStorage() const
    {
        return newMatrixStorage( getNumRows(), getNumColumns() );
    }

    /** Implementation of MatrixStorage::copy for derived class. */

    virtual DenseStorage<ValueType>* copy() const;

    /** Test the storage data for inconsistencies.
     *
     *  @param[in] msg is additional message string that should be used to identify calling routine
     *  @throws Exception in case of any inconsistency.
     */
    void check( const char* msg ) const;

    void transposeImpl();

    /** Getter routine for the data type of the elements. */

    virtual common::ScalarType getValueType() const;

    /** Getter routine for the enum value that stands for this format. */

    virtual Format getFormat() const;

    /**
     * @brief sets data for identity matrix.
     *
     * @param[in] size is the size of the square matrix
     */
    void setIdentity( const IndexType size );

    /**
     *  @brief Sets allocated dense matrix to zero
     */
    void setZero();

    /* ==================================================================== */
    /*  set / get CSR data                                                  */
    /* ==================================================================== */

    /** Implementation of _MatrixStorage::setCSRData for this class.  */

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

    /* Print relevant information about matrix storage format. */

    virtual void writeAt( std::ostream& stream ) const;

    /** Get a value of the matrix.
     *
     * @param[in] i is the row index, 0 <= i < getNumRows
     * @param[in] j is the colum index, 0 <= j < mNumColumns
     *
     * Out-of-range check is enabled for DEBUG version.
     */

    ValueType getValue( const IndexType i, const IndexType j ) const;

    /** Implementation of pure method MatrixStorage<ValueType>::setValue for dense storage */

    void setValue( const IndexType i, const IndexType j, const ValueType val,
                   const common::BinaryOp op = common::BinaryOp::COPY );

    /** Initiate an asynchronous data transfer to a specified location. */

    void prefetch( const hmemo::ContextPtr location ) const;

    /** Will wait for all outstanding asynchronous data transfers. */

    void wait() const;

    /** Getter routine for the matrix data array (write). */

    hmemo::HArray<ValueType>& getData();

    /** Getter routine for the matrix data array (read-only). */

    const hmemo::HArray<ValueType>& getValues() const;

    /** Getter routine for the number of stored values. */

    virtual IndexType getNumValues() const;

    /**
     *  Implementation of pure method of _MatrixStorage::clear
     */
    virtual void clear()
    {
        _MatrixStorage::setDimension( 0, 0 );
        mData.clear();
    }

    /**
     *  IImplementation of pure method of class _MatrixStorage::purge
     */
    virtual void purge()
    {
        _MatrixStorage::setDimension( 0, 0 );
        mData.purge();
    }

    /* ========================================================================= */
    /*       Filling dense storage with assembled COO data                      */
    /* ========================================================================= */

    /**
     *  @brief Override MatrixStorage<ValueType>::fillCOO 
     *
     *  Default implementation (conversion to CSR and back) is very inefficient
     */
    virtual void fillCOO(
        hmemo::HArray<IndexType> ia,
        hmemo::HArray<IndexType> ja,
        hmemo::HArray<ValueType> values,
        const common::BinaryOp op = common::BinaryOp::COPY );

    /******************************************************************
     *   invert                                                        *
     ******************************************************************/

    /** @brief Implementation of MatrixStorage::invert for dense storage
     *
     *  Note: other matrix storage can be aliased with this storage.
     */
    virtual void invert( const MatrixStorage<ValueType>& other );

    /** Implementation of MatrixStorage::matrixTimesVector for Dense */

    virtual void matrixTimesVector(
        hmemo::HArray<ValueType>& result,
        const ValueType alpha,
        const hmemo::HArray<ValueType>& x,
        const ValueType beta,
        const hmemo::HArray<ValueType>& y,
        const common::MatrixOp op ) const;

    /** Implementation of MatrixStorage::matrixTimesTimes for Dense */

    virtual void matrixTimesMatrix(
        const ValueType alpha,
        const MatrixStorage<ValueType>& a,
        const MatrixStorage<ValueType>& b,
        const ValueType beta,
        const MatrixStorage<ValueType>& c );

    /** Implementation of MatrixStorage::jacobiIterate for Dense */

    virtual void jacobiIterate(
        hmemo::HArray<ValueType>& solution,
        const hmemo::HArray<ValueType>& oldSolution,
        const hmemo::HArray<ValueType>& rhs,
        const ValueType omega ) const;

    /** Implementation of MatrixStorage::jacobiIterateHalo for Dense */

    virtual void jacobiIterateHalo(
        hmemo::HArray<ValueType>& localSolution,
        const hmemo::HArray<ValueType>& localDiagonal,
        const hmemo::HArray<ValueType>& haloOldSolution,
        const ValueType omega ) const;

    /** Implementation of element-wise binary operation for dense storage */

    virtual void binaryOp(
        const MatrixStorage<ValueType>& a,
        const common::BinaryOp op,
        const MatrixStorage<ValueType>& b );

    /** Implementation for MatrixStorage::l1Norm */

    virtual RealType<ValueType> l1Norm() const;

    /** Implementation for MatrixStorage::l2Norm */

    virtual RealType<ValueType> l2Norm() const;

    /** Implementation for MatrixStorage::maxNorm */

    virtual RealType<ValueType> maxNorm() const;

    /** Implementation for MatrixStorage::maxDiffNorm */

    virtual RealType<ValueType> maxDiffNorm( const MatrixStorage<ValueType>& other ) const;

    /** Implementation of maxDiffNorm for dense matrices */

    virtual RealType<ValueType> maxDiffNormImpl( const DenseStorage<ValueType>& other ) const;

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

    /** Implementation of pure method MatrixStorage<ValueType>::scaleRows */

    void scaleRows( const hmemo::HArray<ValueType>& values );

    /** Implementation of pure method MatrixStorage<ValueType>::scaleColumns */

    void scaleColumns( const hmemo::HArray<ValueType>& values );

    /** Implementation of pure method.  */

    virtual void scale( const ValueType value );

    /** Implementation of pure method.  */

    virtual void conj();

    virtual size_t getMemoryUsageImpl() const;

    void swap( DenseStorage<ValueType>& other );

    /** Implementation of pure method _MatrixStorage::allocate */

    void allocate( const IndexType numRows, const IndexType numColumns );

    /** Print dense matrix on output; use it carefully */

    virtual void print( std::ostream& ) const;

    using MatrixStorage<ValueType>::prefetch;
    using MatrixStorage<ValueType>::assign;
    using MatrixStorage<ValueType>::getContextPtr;

    using _MatrixStorage::getNumRows;
    using _MatrixStorage::getNumColumns;

protected:

    using MatrixStorage<ValueType>::mContext;

    hmemo::HArray<ValueType> mData;  //!<  matrix data as HArray, stored row-wise

    /** Logger just for this class / matrix format. */

    SCAI_LOG_DECL_STATIC_LOGGER( logger )

private:

    /** Implementation of matrix times matrix for dense matrices. */

    void matrixTimesMatrixDense( const ValueType alpha,
                                 const DenseStorage<ValueType>& a,
                                 const DenseStorage<ValueType>& b,
                                 const ValueType beta,
                                 const DenseStorage<ValueType>& c );

    /** Implementation of elementwise binary operation for dense matrices. */

    void binaryOpDense(
        const DenseStorage<ValueType>& a,
        const common::BinaryOp op,
        const DenseStorage<ValueType>& b );

    /** @brief invert only for DenseStorage. */

    void invertDense( const DenseStorage<ValueType>& other );

    static std::string initTypeName();
};

/* --------------------------------------------------------------------------- */

} /* end namespace lama */

} /* end namespace scai */
