/**
 * @file ELLStorage.hpp
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
 * @brief Definition of a structure for a (non-distributed) ELL sparse matrix.
 * @author Lauretta Schubert
 * @date 25.05.2011
 * @since 1.0.0
 */

#ifndef LAMA_ELLSTORAGE_HPP_
#define LAMA_ELLSTORAGE_HPP_

// for dll_import
#include <lama/config.hpp>

// base classes
#include <lama/storage/CRTPMatrixStorage.hpp>

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
 *  @tparam T is the value type of the matrix values.
 */
template<typename T>
class LAMA_DLL_IMPORTEXPORT ELLStorage: public CRTPMatrixStorage<ELLStorage<T>, T>
{
public:

    typedef T ValueType; //!< This is the type of the matrix values.

    /** get typename of the matrix storage format. */

    static const char* typeName();

    /**
     *  Default constructor, zero matrix.
     */
    ELLStorage();

    /**
     * @brief Creates a sparse matrix with all values set to zero.
     *
     * @param[in] numRows       the number of rows of the matrix
     * @param[in] numColumns    the number of columns of the matrix
     * @param[in] con           the context type where conversion should be done
     */
    ELLStorage( const IndexType numRows, const IndexType numColumns, const ContextType con = Context::Host );

    /** Constructor for ELL storage by corresponding arrays.
     *
     * @param[in] numRows           number of rows
     * @param[in] numColumns        number of columns
     * @param[in] numValuesPerRows  the number of stored elements in the matrix
     * @param[in] ia                row pointer of the input matrix
     * @param[in] ja                column indexes of the input matrix
     * @param[in] values            the data values of the input matrix
     */
    ELLStorage(
        const IndexType numRows,
        const IndexType numColumns,
        const IndexType numValuesPerRows,
        const LAMAArray<IndexType>& ia,
        const LAMAArray<IndexType>& ja,
        const LAMAArray<ValueType>& values );

    /** Override the default copy constructor */

    ELLStorage( const ELLStorage<ValueType>& other );

    /** Copy constructor can take any matrix storage. */

    explicit ELLStorage( const _MatrixStorage& other )
    {
        assign( other );
    }

    /** Copy constructor can take any matrix storage or context. */

    ELLStorage( const _MatrixStorage& other, const ContextPtr context )
    {
        setContext( context );
        assign( other );
    }

    /** Destructor of ELL sparse matrix. */

    virtual ~ELLStorage();

    /** Override the default assignment operator */

    ELLStorage<ValueType>& operator=( const ELLStorage<ValueType>& other );

    /** The assignment operator can be used with any matrix storage.
     *
     *  @param[in] other is the matrix storage that will be assigned
     *
     *  The assignment operator will also carry out implicit type and
     *  and format conversions.
     */

    ELLStorage<ValueType>& operator=( const _MatrixStorage& other );

    /** Implementation of MatrixStorage::copy for derived class. */

    virtual ELLStorage* copy() const;

    /** Implementation of MatrixStorage::create for derived class. */

    virtual ELLStorage* create() const;

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

    virtual MatrixStorageFormat getFormat() const;

    /** Implementation of pure method of class MatrixStorage. */

    virtual void purge();

    /** Implementation of pure method. */

    virtual void setIdentity( const IndexType size );

    /** General routine to build any kind of CSR storage.
     *
     * @param[out] ia is the CSR offset array
     * @param[out] ja is the array with the column indexes (optional)
     * @param[out] values is the array with the non-zero matrix values (optional)
     * @param[in]  loc is the Context where conversion should be done
     */
    template<typename OtherValueType>
    void buildCSR(
        LAMAArray<IndexType>& ia,
        LAMAArray<IndexType>* ja,
        LAMAArray<OtherValueType>* values,
        const ContextPtr loc ) const;

    /**
     * @brief fills ELLPACK sparse matrix by csr sparse data.
     *
     * @param[in] numRows    number of rows
     * @param[in] numColumns number of columns
     * @param[in] numValues  the number of stored elements in the matrix
     * @param[in] ia         row pointer of the input csr sparse matrix
     * @param[in] ja         column indexes of the input csr sparse matrix
     * @param[in] values     the data values of the input csr sparse matrix
     * @param[in]  loc is the Context where conversion should be done
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

    void setELLData(
        const IndexType numRows,
        const IndexType numColumns,
        const IndexType numValuesPerRow,
        const LAMAArray<IndexType>& ia,
        const LAMAArray<IndexType>& ja,
        const _LAMAArray& values );

    /* Print relevant information about matrix storage format. */

    virtual void writeAt( std::ostream& stream ) const;

    /** getter for member variables IA, JA, Data, only const reference */

    const LAMAArray<IndexType>& getIA() const;

    const LAMAArray<IndexType>& getJA() const;

    const LAMAArray<ValueType>& getValues() const;

    /** Getter routine for the number of stored values*/

    virtual IndexType getNumValues() const;

    /** Getter routine for the number of num values per row*/

    IndexType getNumValuesPerRow() const;

    /** Template method for getting row. */

    template<typename OtherType>
    void getRowImpl( LAMAArray<OtherType>& row, const IndexType i ) const;

    /** This method returns the diagonal
     *
     * @param[in] diagonal  is the destination array
     *
     * Calculations are dependent to the diagonal property
     */
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

    /** This method replaces the diagonal by a diagonal value
     *
     * @param[in] scalar  is the source value
     *
     * Calculations are dependent to the diagonal property
     */
    void setDiagonalImpl( const Scalar scalar );

    /******************************************************************
     *  Scaling of elements in a matrix                                *
     ******************************************************************/

    /** Template version used for virtual routine scale with known value type. */

    template<typename OtherType>
    void scaleImpl( const LAMAArray<OtherType>& values );

    /** Implementation of pure method.  */

    void scaleImpl( const Scalar value );

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

    /** Implementation of MatrixStorage::matrixTimesVector for ELL */

    virtual void matrixTimesVector(
        LAMAArrayView<ValueType> result,
        const ValueType alpha,
        const LAMAArrayConstView<ValueType> x,
        const ValueType beta,
        const LAMAArrayConstView<ValueType> y ) const;

    /** Implementation of MatrixStorage::matrixTimesVectorAsync for ELL */

    virtual SyncToken* matrixTimesVectorAsync(
        LAMAArrayView<ValueType> result,
        const ValueType alpha,
        const LAMAArrayConstView<ValueType> x,
        const ValueType beta,
        const LAMAArrayConstView<ValueType> y ) const;

    /** Implementation of MatrixStorage::jacobiIterate for ELL */

    virtual void jacobiIterate(
        LAMAArrayView<ValueType> solution,
        const LAMAArrayConstView<ValueType> oldSolution,
        const LAMAArrayConstView<ValueType> rhs,
        const ValueType omega ) const;

    /** Implementation of MatrixStorage::jacobiIterateHalo for ELL */

    virtual void jacobiIterateHalo(
        LAMAArrayView<ValueType> localSolution,
        const MatrixStorage<ValueType>& localStorage,
        const LAMAArrayConstView<ValueType> haloOldSolution,
        const ValueType omega ) const;

    /** Implementation of MatrixStorage::jacobiIterateHalo for ELL 
     *  @since 1.1.0
     */

    virtual void jacobiIterateHalo(
        LAMAArrayView<ValueType> localSolution,
        const LAMAArray<ValueType>* localDiagonal,
        const LAMAArrayConstView<ValueType> haloOldSolution,
        const ValueType omega ) const;

    /** Implementation of MatrixStorage::matrixTimesMatrix for ELL */

    virtual void matrixTimesMatrix(
        const ValueType alpha,
        const MatrixStorage<ValueType>& a,
        const MatrixStorage<ValueType>& b,
        const ValueType beta,
        const MatrixStorage<ValueType>& c );

    /** Implementation for MatrixStorage::maxNorm */

    virtual ValueType maxNorm() const;

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

    /** Swap this ELL storage data with another ELL storage.
     *
     *  @param[in,out] other is the ELL storage with with this storage is swapped
     *
     *  Note: swap is only possible for two storages of the same format and same type.
     */
    void swap( ELLStorage<ValueType>& other );

    virtual size_t getMemoryUsageImpl() const;

    using MatrixStorage<ValueType>::mNumRows;
    using MatrixStorage<ValueType>::mNumColumns;
    using MatrixStorage<ValueType>::mDiagonalProperty;
    using MatrixStorage<ValueType>::mRowIndexes;
    using MatrixStorage<ValueType>::mCompressThreshold;

    using MatrixStorage<ValueType>::assign;
    using MatrixStorage<ValueType>::prefetch;
    using MatrixStorage<ValueType>::getContextPtr;
    using MatrixStorage<ValueType>::setContext;

    /** Print out the ELL storage on std::out, use only for debug. */

    void print() const;

private:

    IndexType mNumValuesPerRow; //!< number of values in each row

    LAMAArray<IndexType> mIA; //!< size is numRows
    LAMAArray<IndexType> mJA; //!< size is numRows x numValuesPerRow
    LAMAArray<ValueType> mValues; //!< size is numRows x numValuesPerRow

    /** Addressing function for the arrays ia and ja: column-wise */

    inline IndexType ellindex( const IndexType i, const IndexType jj ) const
    {
        return jj * mNumRows + i;
    }

    // Alternative addressing row-wise
    // inline IndexType ellindex( const IndexType i, const IndexType jj ) const { return i * mNumValuesPerRow + jj; }

    /**
     * @brief checks storage data if diagonal property is given
     *
     */
    virtual bool checkDiagonalProperty() const;

    /** Help routine that computes array with row indexes for non-empty rows.
     *  The array is only built if number of non-zero rows is smaller than
     *  a certain percentage ( mThreshold ).
     */

    void buildRowIndexes( const ContextPtr loc );

    LAMA_LOG_DECL_STATIC_LOGGER( logger ) //!< logger for this matrix format

    /** Matrix times Matrix for ELL only */
    void matrixTimesMatrixELL( const ValueType alpha,
                               const ELLStorage<ValueType>& a,
                               const ELLStorage<ValueType>& b );

    /** Matrix + Matrix for ELL only */
    void matrixAddMatrixELL( const ValueType alpha,
                             const ELLStorage<ValueType>& a,
                             const ValueType beta,
                             const ELLStorage<ValueType>& b );
};

} // namespace lama

#endif // LAMA_ELLSTORAGE_HPP_
