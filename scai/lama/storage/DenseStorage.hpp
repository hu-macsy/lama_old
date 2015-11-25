/**
 * @file DenseStorage.hpp
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
 * @brief Definition of a structure for a dense matrix.
 * @author Michael Drost
 * @date 27.04.2011
 * @since 1.0.0
 */

#pragma once

// for dll_import
#include <scai/common/config.hpp>

// base classes
#include <scai/lama/storage/CRTPMatrixStorage.hpp>

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
class COMMON_DLL_IMPORTEXPORT DenseStorageView: public CRTPMatrixStorage<DenseStorageView<ValueType>,ValueType>
{
public:

    typedef ValueType StorageValueType;

    /** get typename of the matrix storage format. */

    static const char* typeName();

    /**
     * @brief Creates a dense matrix with referenced data array
     *
     * @param[in] data            is an array whose values will be copied
     * @param[in] numRows         the number of rows of the matrix
     * @param[in] numColumns      the number of columns of the matrix
     * @param[in] initializedData should be true for existing data to verify good sizes
     */
    DenseStorageView(
        hmemo::HArray<ValueType>& data,
        const IndexType numRows,
        const IndexType numColumns,
        bool initializedData = true );

    /** Destructor of Dense matrix reference. */

    virtual ~DenseStorageView();

    /** Implementation of MatrixStorage::copy for derived class. */

    virtual DenseStorageView<ValueType>* copy() const;

    /** Implementation of MatrixStorage::create for derived class. */

    virtual DenseStorageView<ValueType>* clone() const;

    /** Test the storage data for inconsistencies.
     *
     *  @param[in] msg is additional message string that should be used to identify calling routine
     *  @throws Exception in case of any inconsistency.
     */
    void check( const char* msg ) const;

    /** Override default assignment for more efficient solutions */

    virtual void assign( const _MatrixStorage& other );

    void transposeImpl();
    
    /** Getter routine for the data type of the elements. */

    virtual common::scalar::ScalarType getValueType() const;

    /** Getter routine for the enum value that stands for this format. */

    virtual MatrixStorageFormat getFormat() const;

    /**
     * @brief sets data for identity matrix.
     *
     * @param[in] size is the size of the square matrix
     */
    void setIdentity( const IndexType size );

    /**
     *  @brief Sets current dense matrix to identity
     */
    void setIdentity();

    /**
     *  @brief Sets allocated dense matrix to zero
     */
    void setZero();

    /**
     * @brief fills dense storage by csr sparse data.
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

    /* Print relevant information about matrix storage format. */

    virtual void writeAt( std::ostream& stream ) const;

    /** Get a value of the matrix.
     *
     * @param[in] i is the row index, 0 <= i < mNumRows
     * @param[in] j is the colum index, 0 <= j < mNumRows
     *
     * Out-of-range check is enabled for DEBUG version.
     */

    ValueType getValue( const IndexType i, const IndexType j ) const;

    /** Initiate an asynchronous data transfer to a specified location. */

    void prefetch( const hmemo::ContextPtr location ) const;

    /** Will wait for all outstanding asynchronous data transfers. */

    void wait() const;

    /** Getter routine for the matrix data array (write). */

    hmemo::HArray<ValueType>& getData();

    /** Getter routine for the matrix data array (read-only). */

    const hmemo::HArray<ValueType>& getData() const;

    /** Getter routine for the number of stored values. */

    virtual IndexType getNumValues() const;

    /**
     *  Implementation of pure method of _MatrixStorage::clear
     */
    virtual void clear()
    {
        mNumRows = 0;
        mNumColumns = 0;
        mData.clear();
        mDiagonalProperty = checkDiagonalProperty();
    }

    /**
     *  IImplementation of pure method of class _MatrixStorage::purge
     */
    virtual void purge()
    {
        mNumRows = 0;
        mNumColumns = 0;
        mData.purge();
        mDiagonalProperty = checkDiagonalProperty();
    }

    /** Help routine to build any kind of CSR storage.
     *
     *
     *  @param[out] csrIA is the CSR offset array
     *  @param[out] csrJA is the array with the column indexes (optional)
     *  @param[out] csrValues is the array with the non-zero matrix values (optional)
     *  @param[in]  loc is the Context where conversion should be done
     */
    template<typename OtherValueType>
    void buildCSR(
        hmemo::HArray<IndexType>& csrIA,
        hmemo::HArray<IndexType>* csrJA,
        hmemo::HArray<OtherValueType>* csrValues,
        const hmemo::ContextPtr loc ) const;

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
        const hmemo::HArray<ValueType>& y ) const;

    /** Implementation of MatrixStorage::vectorTimesMatrix for Dense */
    /** since 1.0.1 */

    virtual void vectorTimesMatrix(
        hmemo::HArray<ValueType>& result,
        const ValueType alpha,
        const hmemo::HArray<ValueType>& x,
        const ValueType beta,
        const hmemo::HArray<ValueType>& y ) const;

    /** Implementation of MatrixStorage::matrixTimesTimes for Dense */

    virtual void matrixTimesMatrix(
        const ValueType alpha,
        const MatrixStorage<ValueType>& a,
        const MatrixStorage<ValueType>& b,
        const ValueType beta,
        const MatrixStorage<ValueType>& c );

    /** Implementation for MatrixStorage::l1Norm */

    virtual ValueType l1Norm() const;

    /** Implementation for MatrixStorage::l2Norm */

    virtual ValueType l2Norm() const;

    /** Implementation for MatrixStorage::maxNorm */

    virtual ValueType maxNorm() const;

    /** Implementation for MatrixStorage::maxDiffNorm */

    virtual ValueType maxDiffNorm( const MatrixStorage<ValueType>& other ) const;

    /** Implementation of maxDiffNorm for dense matrices */

    virtual ValueType maxDiffNormImpl( const DenseStorageView<ValueType>& other ) const;

    /** Template method for getting row. */

    template<typename OtherType>
    void getRowImpl( hmemo::HArray<OtherType>& row, const IndexType i ) const;

    /** Template method for getting diagonal. */

    template<typename OtherType>
    void getDiagonalImpl( hmemo::HArray<OtherType>& diagonal ) const;

    /** This method replaces the diagonal
     *
     * @param[in] diagonal  is the source array
     *
     * Calculations are dependent to the diagonal property
     */
    template<typename OtherType>
    void setDiagonalImpl( const hmemo::HArray<OtherType>& diagonal );

    /** Implementation of pure method. */

    void setDiagonalImpl( const ValueType value );

    /******************************************************************
     *  Scaling of elements in a matrix                                *
     ******************************************************************/

    /** Template version used for virtual routine scale with known value type. */

    template<typename OtherType>
    void scaleImpl( const hmemo::HArray<OtherType>& values );

    /** Implementation of pure method.  */

    virtual void scaleImpl( const ValueType value );

    virtual size_t getMemoryUsageImpl() const;

    void swap( DenseStorageView<ValueType>& other );

    void allocate( const IndexType numRows, const IndexType numColumns );

    /** Print dense matrix on output; use it carefully */

    void print() const;

    using MatrixStorage<ValueType>::prefetch;
    using MatrixStorage<ValueType>::assign;
    using MatrixStorage<ValueType>::getTypeName;

protected:

    using MatrixStorage<ValueType>::mNumRows;
    using MatrixStorage<ValueType>::mNumColumns;
    using MatrixStorage<ValueType>::mDiagonalProperty;
    using MatrixStorage<ValueType>::mContext;

    hmemo::HArray<ValueType>& mData; //!<  Reference to the matrix value array

    /** Logger just for this class / matrix format. */

    SCAI_LOG_DECL_STATIC_LOGGER( logger )

    /** Override MatrixStorage::checkDiagonalProperty method. */

    virtual    bool checkDiagonalProperty() const;

private:

    /** Disable default constructor. */

    DenseStorageView();

    template<typename OtherType>
    void assignDenseStorageImpl( const DenseStorageView<OtherType>& otherDenseStorage );

    /** Implementation of matrix times matrix for dense matrices. */

    void matrixTimesMatrixDense( const ValueType alpha,
                                 const DenseStorageView<ValueType>& a,
                                 const DenseStorageView<ValueType>& b,
                                 const ValueType beta,
                                 const DenseStorageView<ValueType>& c );

    /** @brief invert only for DenseStorage. */

    void invertDense( const DenseStorageView<ValueType>& other );
};

/* --------------------------------------------------------------------------- */

/** DenseStorage is a DenseStorageView that manages its own data.
 *
 *  @tparam ValueType is the value type of the matrix values.
 */
template<typename ValueType>
class COMMON_DLL_IMPORTEXPORT DenseStorage: public DenseStorageView<ValueType>
{
public:

    /** @brief Getter of type name of the matrix storage format. */

    static const char* typeName();

    /** @brief Default constructor for a dense storage of 0 x 0  */

    DenseStorage();

    /**
     * @brief Creates a dense matrix with all values set to zero.
     *
     * @param[in] numRows       the number of rows of the matrix
     * @param[in] numColumns    the number of columns of the matrix
     */
    DenseStorage( const IndexType numRows, const IndexType numColumns );

    /**
     * @brief Creates a dense matrix with all values copied
     *
     * @param[in] data          is an array whose values will be copied
     * @param[in] numRows       the number of rows of the matrix
     * @param[in] numColumns    the number of columns of the matrix
     */
    DenseStorage( const hmemo::HArray<ValueType>& data, const IndexType numRows, const IndexType numColumns );

    /** Default copy constructor is overridden */

    DenseStorage( const DenseStorage<ValueType>& other );

    /** Copy constructor that handles also type and format conversion. */

    explicit DenseStorage( const _MatrixStorage& other );

    /** Copy constructor can take any matrix storage or context. */

    DenseStorage( const _MatrixStorage& other, const hmemo::ContextPtr context );

    /** Override MatrixStorage::getTypeName() */

    const char* getTypeName() const;

    DenseStorage<ValueType>& operator=( const _MatrixStorage& other );

    DenseStorage<ValueType>& operator=( const DenseStorage<ValueType>& other );

    /** Destructor of Dense matrix. */

    virtual ~DenseStorage();

    using MatrixStorage<ValueType>::mNumRows;
    using MatrixStorage<ValueType>::mNumColumns;
    using MatrixStorage<ValueType>::assign;

    using DenseStorageView<ValueType>::mData;

private:

    hmemo::HArray<ValueType> mDataArray; //!<  matrix values, size is mNumRows x mNumColumns

};

} /* end namespace lama */

} /* end namespace scai */
