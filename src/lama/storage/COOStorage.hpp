/**
 * @file COOStorage.hpp
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
 * @brief Definition of a matrix storage class for a (non-distributed) sparse matrix
 *        represented by coordinate lists.
 * @author Paul Börding, Thomas Brandes
 * @date 11.06.2011
 * @since 1.0.0
 */
#ifndef LAMA_COOSTORAGE_HPP_
#define LAMA_COOSTORAGE_HPP_

// for dll_import
#include <lama/config.hpp>

// base classes
#include <lama/storage/CRTPMatrixStorage.hpp>

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
 * @tparam T is the value type of the matrix values.
 *
 * Note: Copy constructor and operator= are supported by default implementations,
 *       copies of LAMAArray will be deep copies.
 */

template<typename T>
class LAMA_DLL_IMPORTEXPORT COOStorage: public CRTPMatrixStorage<COOStorage<T>,T>
{
public:

    typedef T ValueType; //!< This is the type of the matrix values.

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
        const LAMAArray<IndexType>& ia,
        const LAMAArray<IndexType>& ja,
        const LAMAArray<ValueType>& values );

    /** Default constructor, same as COOStorage(0, 0). */

    COOStorage();

    /** Default copy constructor is overridden */

    COOStorage( const COOStorage<ValueType>& other );

    /** Copy constructor that handles also type and format conversion. */

    explicit COOStorage( const _MatrixStorage& other )
    {
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

    virtual MatrixStorageFormat getFormat() const;

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

    /** Implementation of pure method for COO matrix. */

    virtual void setIdentity( const IndexType size );

    /** Help routine to build any kind of CSR storage. */

    template<typename OtherValueType>
    void buildCSR(
        LAMAArray<IndexType>& ia,
        LAMAArray<IndexType>* ja,
        LAMAArray<OtherValueType>* values,
        const ContextPtr /* loc */) const;

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
        const LAMAArray<IndexType>& ia,
        const LAMAArray<IndexType>& ja,
        const _LAMAArray& values );

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
        const LAMAArray<IndexType>& ia,
        const LAMAArray<IndexType>& ja,
        const LAMAArray<OtherValueType>& values,
        const ContextPtr loc ) __attribute__( ( noinline ) );

    /* Print relevant information about matrix storage format. */

    virtual void writeAt( std::ostream& stream ) const;

    /** Getter for member variable IA, only const reference. */

    const LAMAArray<IndexType>& getIA() const;

    /** Getter for member variable JA, only const reference. */

    const LAMAArray<IndexType>& getJA() const;

    /** Getter for member variable values, only const reference. */

    const LAMAArray<ValueType>& getValues() const;

    /** Getter for the number of stored values. */

    virtual IndexType getNumValues() const;

    /** Template method for getting row. */

    template<typename OtherType>
    void getRowImpl( LAMAArray<OtherType>& row, const IndexType i ) const __attribute( ( noinline ) );

    /** This method returns the diagonal
     *
     * @param[in] diagonal  is the destination array
     *
     * Calculations are dependent to the diagonal property
     */
    template<typename OtherType>
    void getDiagonalImpl( LAMAArray<OtherType>& diagonal ) const __attribute( ( noinline ) );

    /** Template version used for virtual routine setDiagonal with known value type. */

    template<typename OtherType>
    void setDiagonalImpl( const LAMAArray<OtherType>& diagonal ) __attribute( ( noinline ) );

    /** Implementation of pure method. */

    virtual void setDiagonalImpl( const Scalar scalar );

    /** Implementation of MatrixStorage::copy for derived class. */

    virtual COOStorage* copy() const;

    /** Implementation of MatrixStorage::create for derived class. */

    virtual COOStorage* create() const;

    /******************************************************************
     *  Matrix times Vector                                            *
     ******************************************************************/

    /** Implementation of MatrixStorage::matrixTimesVector for COO */

    virtual void matrixTimesVector(
        LAMAArrayView<ValueType> result,
        const ValueType alpha,
        const LAMAArrayConstView<ValueType> x,
        const ValueType beta,
        const LAMAArrayConstView<ValueType> y ) const;

    /** Implementation of MatrixStorage::vectorTimesMatrix for COO */
    /** since 1.0.1 */

    virtual void vectorTimesMatrix(
        LAMAArray<ValueType>& result,
        const ValueType alpha,
        const LAMAArray<ValueType>& x,
        const ValueType beta,
        const LAMAArray<ValueType>& y ) const;

    /** Implementation of MatrixStorage::matrixTimesVectorAsync for COO */

    virtual std::auto_ptr<SyncToken> matrixTimesVectorAsyncToDo(
        LAMAArrayView<ValueType> result,
        const ValueType alpha,
        const LAMAArrayConstView<ValueType> x,
        const ValueType beta,
        const LAMAArrayConstView<ValueType> y ) const;

    /** Implementation of MatrixStorage::vectorTimesMatrixAsync for CSR */
    /** since 1.0.1 */

    virtual SyncToken* vectorTimesMatrixAsync(
        LAMAArray<ValueType>& result,
        const ValueType alpha,
        const LAMAArray<ValueType>& x,
        const ValueType beta,
        const LAMAArray<ValueType>& y ) const;

    /** solution = xxx */

    virtual void jacobiIterate(
        LAMAArray<ValueType>& solution,
        const LAMAArray<ValueType>& oldSolution,
        const LAMAArray<ValueType>& rhs,
        const ValueType omega ) const;

    /******************************************************************
     *  Scaling of elements in a matrix                                *
     ******************************************************************/

    /** Template version used for virtual routine scale with known value type. */

    template<typename OtherType>
    void scaleImpl( const LAMAArray<OtherType>& values ) __attribute( ( noinline ) );

    /** Implementation of pure method.  */

    void scaleImpl( const Scalar value );

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

    /** Initiate an asynchronous data transfer to a specified location. */

    void prefetch( const ContextPtr location ) const;

    /** Will wait for all outstanding asynchronous data transfers. */

    void wait() const;

    /** Swaps this with other.
     * @param[in,out] other the COOStorage to swap this with
     */
    void swap( COOStorage<ValueType>& other );

    virtual size_t getMemoryUsageImpl() const;

    using MatrixStorage<ValueType>::mNumRows;
    using MatrixStorage<ValueType>::mNumColumns;
    using MatrixStorage<ValueType>::mDiagonalProperty;
    using MatrixStorage<ValueType>::mRowIndexes;
    using MatrixStorage<ValueType>::mCompressThreshold;

    using MatrixStorage<ValueType>::assign;
    using MatrixStorage<ValueType>::prefetch;
    using MatrixStorage<ValueType>::getContextPtr;

    IndexType mNumValues; //!< number of non-zero values (+ optionally zeros in diagonal)

    LAMAArray<IndexType> mIA; //!< row indices, size is mNumValues
    LAMAArray<IndexType> mJA; //!< column indices, size is mNumValues
    LAMAArray<ValueType> mValues; //!< non-zero values (+ optionally zeros in diagonal), size is mNumValues

private:

    /** Function that checks the diagonal property of the COO matrix. If on
     *  the first column index for each row is the row (diagonal element).
     */
    virtual bool checkDiagonalProperty() const;

    // values are stored in row-major order

    // TODO: index at irow / icolumn, not supported yet

    // inline IndexType index(IndexType irow, IndexType icolumn) const { return icolumn * mNumRows + irow; }

    LAMA_LOG_DECL_STATIC_LOGGER( logger )
};

} // namespace lama

#endif // LAMA_COOSTORAGE_HPP_
