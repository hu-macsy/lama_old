/**
 * @file DenseStorage.hpp
 *
 * @license
 * Copyright (c) 2009-2013
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
#ifndef LAMA_DENSESTORAGE_HPP_
#define LAMA_DENSESTORAGE_HPP_

// for dll_import
#include <lama/config.hpp>

// base classes
#include <lama/storage/CRTPMatrixStorage.hpp>

// others
#include <lama/LAMAArray.hpp>
#include <lama/HostReadAccess.hpp>
#include <lama/HostWriteAccess.hpp>
#include <lama/Scalar.hpp>

namespace lama
{

template<typename T> class DenseStorage;

/**
 * Storage format for a Dense matrix.
 *
 *  @tparam T is the value type of the matrix values.
 */
template<typename T>
class LAMA_DLL_IMPORTEXPORT DenseStorageView: public CRTPMatrixStorage<DenseStorageView<T>,T>
{
public:

    typedef T ValueType; //!< This is the type of the matrix values.

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
        LAMAArray<ValueType>& data,
        const IndexType numRows,
        const IndexType numColumns,
        bool initializedData = true );

    /** Destructor of Dense matrix reference. */

    virtual ~DenseStorageView();

    /** Implementation of MatrixStorage::copy for derived class. */

    virtual DenseStorageView<ValueType>* copy() const;

    /** Implementation of MatrixStorage::create for derived class. */

    virtual DenseStorageView<ValueType>* create() const;

    /** Test the storage data for inconsistencies.
     *
     *  @param[in] msg is additional message string that should be used to identify calling routine
     *  @throws Exception in case of any inconsistency.
     */
    void check( const char* msg ) const;

    /** Override default assignment for more efficient solutions */

    virtual void assign( const _MatrixStorage& other );

    /** Getter routine for the data type of the elements. */

    virtual Scalar::ScalarType getValueType() const;

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
        const LAMAArray<IndexType>& ia,
        const LAMAArray<IndexType>& ja,
        const LAMAArray<OtherValueType>& values,
        const ContextPtr loc );

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

    void prefetch( const ContextPtr location ) const;

    /** Will wait for all outstanding asynchronous data transfers. */

    void wait() const;

    /** Getter routine for the matrix data array (write). */

    LAMAArray<ValueType>& getData();

    /** Getter routine for the matrix data array (read-only). */

    const LAMAArray<ValueType>& getData() const;

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
     *  @param[out] ia is the CSR offset array
     *  @param[out] ja is the array with the column indexes (optional)
     *  @param[out] values is the array with the non-zero matrix values (optional)
     *  @param[in]  loc is the Context where conversion should be done
     */
    template<typename OtherValueType>
    void buildCSR(
        LAMAArray<IndexType>& ia,
        LAMAArray<IndexType>* ja,
        LAMAArray<OtherValueType>* values,
        const ContextPtr loc ) const;

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
        LAMAArrayView<ValueType> result,
        const ValueType alpha,
        const LAMAArrayConstView<ValueType> x,
        const ValueType beta,
        const LAMAArrayConstView<ValueType> y ) const;

    /** Implementation of MatrixStorage::vectorTimesMatrix for Dense */
    /** since 1.0.1 */

    virtual void vectorTimesMatrix(
        LAMAArray<ValueType>& result,
        const ValueType alpha,
        const LAMAArray<ValueType>& x,
        const ValueType beta,
        const LAMAArray<ValueType>& y ) const;

    /** Implementation of MatrixStorage::matrixTimesTimes for Dense */

    virtual void matrixTimesMatrix(
        const ValueType alpha,
        const MatrixStorage<ValueType>& a,
        const MatrixStorage<ValueType>& b,
        const ValueType beta,
        const MatrixStorage<ValueType>& c );

    /** Implementation for MatrixStorage::maxNorm */

    virtual ValueType maxNorm() const;

    /** Implementation for MatrixStorage::maxDiffNorm */

    virtual ValueType maxDiffNorm( const MatrixStorage<ValueType>& other ) const;

    /** Implementation of maxDiffNorm for dense matrices */

    virtual ValueType maxDiffNormImpl( const DenseStorageView<ValueType>& other ) const;

    /** Template method for getting row. */

    template<typename OtherType>
    void getRowImpl( LAMAArray<OtherType>& row, const IndexType i ) const;

    /** Template method for getting diagonal. */

    template<typename OtherType>
    void getDiagonalImpl( LAMAArray<OtherType>& diagonal ) const;

    /** This method replaces the diagonal
     *
     * @param[in] diagonal  is the source array
     *
     * Calculations are dependent to the diagonal property
     */
    template<typename OtherType>
    void setDiagonalImpl( const LAMAArray<OtherType>& diagonal );

    /** Implementation of pure method. */

    void setDiagonalImpl( const Scalar scalar );

    /******************************************************************
     *  Scaling of elements in a matrix                                *
     ******************************************************************/

    /** Template version used for virtual routine scale with known value type. */

    template<typename OtherType>
    void scaleImpl( const LAMAArray<OtherType>& values );

    /** Implementation of pure method.  */

    virtual void scaleImpl( const Scalar value );

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

    LAMAArray<ValueType>& mData; //!<  Reference to the matrix value array

    /** Logger just for this class / matrix format. */

    LAMA_LOG_DECL_STATIC_LOGGER( logger )

    /** Override MatrixStorage::checkDiagonalProperty method. */

    virtual bool checkDiagonalProperty() const;

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
 *  @tparam T is the value type of the matrix values.
 */
template<typename T>
class LAMA_DLL_IMPORTEXPORT DenseStorage: public DenseStorageView<T>
{
public:

    typedef T ValueType; //!< This is the type of the matrix values.

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
    DenseStorage( const LAMAArray<ValueType>& data, const IndexType numRows, const IndexType numColumns );

    /** Default copy constructor is overridden */

    DenseStorage( const DenseStorage<ValueType>& other );

    /** Copy constructor that handles also type and format conversion. */

    explicit DenseStorage( const _MatrixStorage& other );

    /** Copy constructor can take any matrix storage or context. */

    DenseStorage( const _MatrixStorage& other, const ContextPtr context );

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

    LAMAArray<ValueType> mDataArray; //!<  matrix values, size is mNumRows x mNumColumns

};

}

#endif // LAMA_DENSESTORAGE_HPP_
