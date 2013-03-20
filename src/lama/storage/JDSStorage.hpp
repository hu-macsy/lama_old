/**
 * @file JDSStorage.hpp
 *
 * @license
 * Copyright (c) 2011
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
 * @brief Definition of a structure for a (non-distributed) JDS sparse matrix.
 * @author Thomas Brandes
 * @date 11.06.2011
 * $Id$
 */
#ifndef LAMA_JDSSTORAGE_HPP_
#define LAMA_JDSSTORAGE_HPP_

// for dll_import
#include <lama/config.hpp>

// base classes
#include <lama/storage/CRTPMatrixStorage.hpp>

// others
#include <lama/LAMAArray.hpp>
#include <lama/Scalar.hpp>

#include <lama/distribution/Distribution.hpp>

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
 *  @tparam T is the value type for the matrix elements.
 *
 *  Note: default copy constructor and assignment operator can be used.
 */

template<typename T>
class LAMA_DLL_IMPORTEXPORT JDSStorage: public CRTPMatrixStorage<JDSStorage<T>,T>
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
    JDSStorage( const IndexType numRows, const IndexType numColumns );

    /** Constructor for JDS storage by corresponding arrays. */

    JDSStorage(
        const IndexType numRows,
        const IndexType numColumns,
        const IndexType numValues,
        const IndexType numDiagonals,
        const LAMAArray<IndexType>& dlg,
        const LAMAArray<IndexType>& ilg,
        const LAMAArray<IndexType>& perm,
        const LAMAArray<IndexType>& ja,
        const LAMAArray<ValueType>& values );

    /** Default constructor, same as JDSStorage(0, 0). */

    JDSStorage();

    /** Default copy constructor is overridden */

    JDSStorage( const JDSStorage<ValueType>& other );

    /** Copy constructor that handles also type and format conversion. */

    explicit JDSStorage( const _MatrixStorage& other );

    /** Copy constructor can take any matrix storage or context. */

    JDSStorage( const _MatrixStorage& other, const ContextPtr context );

    /** Default assignment operator is overridden */

    JDSStorage<ValueType>& operator=( const JDSStorage<ValueType>& other );

    /** Assignment operator with same semantic as assign. */

    JDSStorage<ValueType>& operator=( const _MatrixStorage& other );

    /** Implementation of MatrixStorage::copy for derived class. */

    virtual JDSStorage* copy() const;

    /** Implementation of MatrixStorage::create for derived class. */

    virtual JDSStorage* create() const;

    /** Implementation of pure method for _MatrixStorage. */

    virtual void clear();

    /** Destructor of JDS sparse matrix. */

    virtual ~JDSStorage();

    /** Test the storage data for inconsistencies.
     *
     *  @throw Exception in case of any inconsistency.
     */
    void check( const char* msg ) const;

    /** @brief Allocate arrays ilg and perm
     */
    void allocate( const IndexType numRows, const IndexType numColumns );

    /** Getter routine for the enum value that stands for this format. */

    virtual MatrixStorageFormat getFormat() const;

    /** Implementation of pure method of class MatrixStorage. */

    virtual void purge();

    /** Implementation of pure method for JDS matrix. */

    virtual void setIdentity( const IndexType size );

    /** General routine to build any kind of CSR storage.
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
        const ContextPtr /* loc */) const;

    /**
     * Template version with given value type.
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

    /** Implementation of MatrixStorage::matrixTimesVector for JDS */

    virtual void matrixTimesVector(
        LAMAArrayView<ValueType> result,
        const ValueType alpha,
        const LAMAArrayConstView<ValueType> x,
        const ValueType beta,
        const LAMAArrayConstView<ValueType> y ) const;

    /** Implementation of MatrixStorage::jacobiIterate for JDS */

    virtual void jacobiIterate(
        LAMAArrayView<ValueType> solution,
        const LAMAArrayConstView<ValueType> oldSolution,
        const LAMAArrayConstView<ValueType> rhs,
        const ValueType omega ) const;

    /** Implementation of MatrixStorage::jacobiIterateHalo for JDS */

    virtual void jacobiIterateHalo(
        LAMAArrayView<ValueType> localSolution,
        const MatrixStorage<ValueType>& localStorage,
        const LAMAArrayConstView<ValueType> haloOldSolution,
        const ValueType omega ) const;

    /** Print relevant information about matrix storage format. */

    virtual void writeAt( std::ostream& stream ) const;

    const LAMAArray<IndexType>& getJA() const;

    const LAMAArray<IndexType>& getPerm() const;

    const LAMAArray<IndexType>& getDlg() const;

    const LAMAArray<IndexType>& getIlg() const;

    const LAMAArray<ValueType>& getValues() const;

    IndexType getNumDiagonals() const;

    /** Override _MatrixStorage::getNumValues with a more efficient routine for JDS. */

    virtual IndexType getNumValues() const;

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

    void setDiagonalImpl( const Scalar );

    /******************************************************************
     *  Scaling of elements in a matrix                                *
     ******************************************************************/

    /** Template version used for virtual routine scale with known value type. */

    template<typename OtherType>
    void scaleImpl( const LAMAArray<OtherType>& values );

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

    ValueType getValue( IndexType i, IndexType j ) const;

    /** Initiate an asynchronous data transfer to a specified location. */

    void prefetch( const ContextPtr location ) const;

    /** Will wait for all outstanding asynchronous data transfers. */

    void wait() const;

    /** Swaps this with other.
     * @param[in,out] other the JDSStorage to swap this with
     */
    void swap( JDSStorage<ValueType>& other );

    virtual size_t getMemoryUsageImpl() const;

    using MatrixStorage<ValueType>::assign;
    using MatrixStorage<ValueType>::prefetch;
    using MatrixStorage<ValueType>::getContext;
    using MatrixStorage<ValueType>::getContextPtr;

protected:

    using MatrixStorage<ValueType>::mNumRows;
    using MatrixStorage<ValueType>::mNumColumns;
    using MatrixStorage<ValueType>::mDiagonalProperty;
    using MatrixStorage<ValueType>::mRowIndexes;
    using MatrixStorage<ValueType>::mCompressThreshold;

    IndexType mNumDiagonals; //!< number of jagged diagonals (equals length of the longest row)
    IndexType mNumValues; //!< number of non-zero values (+ optionally zeros in diagonal)

    LAMAArray<IndexType> mDlg; //!< number of values in each column, size is mNumDiagonals
    LAMAArray<IndexType> mIlg; //!< number of values in each row, size is mNumRows
    LAMAArray<IndexType> mPerm; //!< position of each row in original matrix, size is mNumRows
    LAMAArray<IndexType> mJa; //!< column indices, size is mNumValues
    LAMAArray<ValueType> mValues; //!< non-zero values (+ optionally zeros in diagonal), size is mNumValues

private:

    /** Assigment with same storage type. */

    void assignJDS( const JDSStorage<ValueType>& other );

    /**
     *  Mandatory routine to check storage data for the diagonal property.
     */
    virtual bool checkDiagonalProperty() const;

    /**
     *  @brief Help routine for stable sort of ilg and setting up perm correctly
     */
    void sortRows( ContextPtr loc );

    /** Help routine that sets up mDlg and allocates mJa and mValues after mIlg defined. */

    void setupData( ContextPtr loc );

    /**  Copy JDS storage but remove all entries where column index
     *   is equal nIndex.
     */

    /*
     * NOT USED ANYMORE
     void setCompressJDS( const JDSStorage<ValueType>& other );
     */

    void print() const;

    LAMA_LOG_DECL_STATIC_LOGGER( logger ); //!< logger for this matrix format
};

} // namespace lama

#endif // LAMA_JDSSTORAGE_HPP_
