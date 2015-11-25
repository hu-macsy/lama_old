/**
 * @file JDSStorage.hpp
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
 * @brief Definition of a structure for a (non-distributed) JDS sparse matrix.
 * @author Thomas Brandes
 * @date 11.06.2011
 * @since 1.0.0
 */

#pragma once

// for dll_import
#include <scai/common/config.hpp>

// base classes
#include <scai/lama/storage/CRTPMatrixStorage.hpp>

// local scai libraries
#include <scai/hmemo/HArray.hpp>

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
class COMMON_DLL_IMPORTEXPORT JDSStorage: public CRTPMatrixStorage<JDSStorage<ValueType>,ValueType>
{
public:

    typedef ValueType StorageValueType;

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
        const hmemo::HArray<IndexType>& dlg,
        const hmemo::HArray<IndexType>& ilg,
        const hmemo::HArray<IndexType>& perm,
        const hmemo::HArray<IndexType>& ja,
        const hmemo::HArray<ValueType>& values );

    /** Default constructor, same as JDSStorage(0, 0). */

    JDSStorage();

    /** Default copy constructor is overridden */

    JDSStorage( const JDSStorage<ValueType>& other );

    /** Copy constructor that handles also type and format conversion. */

    explicit JDSStorage( const _MatrixStorage& other );

    /** Copy constructor can take any matrix storage or context. */

    JDSStorage( const _MatrixStorage& other, const hmemo::ContextPtr context )
    {
        setContextPtr( context );
        assign( other );
    }

    /** Default assignment operator is overridden */

    JDSStorage<ValueType>& operator=( const JDSStorage<ValueType>& other );

    /** Assignment operator with same semantic as assign. */

    JDSStorage<ValueType>& operator=( const _MatrixStorage& other );

    /** Implementation of MatrixStorage::copy for derived class. */

    virtual JDSStorage* copy() const;

    /** Implementation of MatrixStorage::create for derived class. */

    virtual JDSStorage* clone() const;

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
     *  @param[in]  context is the Context where conversion should be done
     */

    template<typename OtherValueType>
    void buildCSR(
        hmemo::HArray<IndexType>& ia,
        hmemo::HArray<IndexType>* ja,
        hmemo::HArray<OtherValueType>* values,
        const hmemo::ContextPtr context ) const;

    /**
     * Template version with given value type.
     */
    template<typename OtherValueType>
    void setCSRDataImpl(
        const IndexType numRows,
        const IndexType numColumns,
        const IndexType numValues,
        const hmemo::HArray<IndexType>& ia,
        const hmemo::HArray<IndexType>& ja,
        const hmemo::HArray<OtherValueType>& values,
        const hmemo::ContextPtr context );

     /**
      * Fill up a JDS storage with the given arrays.
      */
    void setJDSData(
        const IndexType numRows,
        const IndexType numColumns,
        const IndexType numValues,
        const IndexType numDiagonals,
        const hmemo::HArray<IndexType>& dlg,
        const hmemo::HArray<IndexType>& ilg,
        const hmemo::HArray<IndexType>& perm,
        const hmemo::HArray<IndexType>& ja,
        const hmemo::ContextArray& values );

    /** Implementation of MatrixStorage::matrixTimesVector for JDS */

    virtual void matrixTimesVector(
        hmemo::HArray<ValueType>& result,
        const ValueType alpha,
        const hmemo::HArray<ValueType>& x,
        const ValueType beta,
        const hmemo::HArray<ValueType>& y ) const;

    /** Implementation of MatrixStorage::vectorTimesMatrix for JDS */
    /** since 1.0.1 */

    virtual void vectorTimesMatrix(
        hmemo::HArray<ValueType>& result,
        const ValueType alpha,
        const hmemo::HArray<ValueType>& x,
        const ValueType beta,
        const hmemo::HArray<ValueType>& y ) const;

    /** Implementation of MatrixStorage::matrixTimesVectorAsync for JDS */

    virtual tasking::SyncToken* matrixTimesVectorAsync(
        hmemo::HArray<ValueType>& result,
        const ValueType alpha,
        const hmemo::HArray<ValueType>& x,
        const ValueType beta,
        const hmemo::HArray<ValueType>& y ) const;

    /** Implementation of MatrixStorage::vectorTimesMatrixAsync for JDS */
    /** since 1.0.1 */

    virtual tasking::SyncToken* vectorTimesMatrixAsync(
        hmemo::HArray<ValueType>& result,
        const ValueType alpha,
        const hmemo::HArray<ValueType>& x,
        const ValueType beta,
        const hmemo::HArray<ValueType>& y ) const;

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
        const MatrixStorage<ValueType>& localStorage,
        const hmemo::HArray<ValueType>& haloOldSolution,
        const ValueType omega ) const;

    /** Implementation of MatrixStorage::jacobiIterateHalo for JDS
     *  @since 1.1.0
     */

    virtual void jacobiIterateHalo(
        hmemo::HArray<ValueType>& localSolution,
        const hmemo::HArray<ValueType>& localDiagonal,
        const hmemo::HArray<ValueType>& haloOldSolution,
        const ValueType omega ) const;

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

    /** Template method for getting row. */

    template<typename OtherType>
    void getRowImpl( hmemo::HArray<OtherType>& row, const IndexType i ) const __attribute( ( noinline ) );

    /** This method returns the diagonal
     *
     * @param[in] diagonal  is the destination array
     *
     * Calculations are dependent to the diagonal property
     */
    template<typename OtherType>
    void getDiagonalImpl( hmemo::HArray<OtherType>& diagonal ) const __attribute( ( noinline ) );

    /** This method replaces the diagonal
     *
     * @param[in] diagonal  is the source array
     *
     * Calculations are dependent to the diagonal property
     */
    template<typename OtherType>
    void setDiagonalImpl( const hmemo::HArray<OtherType>& diagonal ) __attribute( ( noinline ) );

    void setDiagonalImpl( const ValueType value );

    /******************************************************************
     *  Scaling of elements in a matrix                                *
     ******************************************************************/

    /** Template version used for virtual routine scale with known value type. */

    template<typename OtherType>
    void scaleImpl( const hmemo::HArray<OtherType>& values ) __attribute( ( noinline ) );

    /** Implementation of pure method.  */

    void scaleImpl( const ValueType value );

    /** Implementation for MatrixStorage::l1Norm */

    virtual ValueType l1Norm() const;

    /** Implementation for MatrixStorage::l2Norm */

    virtual ValueType l2Norm() const;

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

    void prefetch( const hmemo::ContextPtr context ) const;

    /** Will wait for all outstanding asynchronous data transfers. */

    void wait() const;

    /** Swaps this with other.
     * @param[in,out] other the JDSStorage to swap this with
     */
    void swap( JDSStorage<ValueType>& other );

    virtual size_t getMemoryUsageImpl() const;

    using MatrixStorage<ValueType>::assign;
    using MatrixStorage<ValueType>::prefetch;
    using MatrixStorage<ValueType>::getContextPtr;
    using MatrixStorage<ValueType>::setContextPtr;

protected:

    using MatrixStorage<ValueType>::mNumRows;
    using MatrixStorage<ValueType>::mNumColumns;
    using MatrixStorage<ValueType>::mDiagonalProperty;
    using MatrixStorage<ValueType>::mRowIndexes;
    using MatrixStorage<ValueType>::mCompressThreshold;

    IndexType mNumDiagonals; //!< number of jagged diagonals (equals length of the longest row)
    IndexType mNumValues; //!< number of non-zero values (+ optionally zeros in diagonal)

    hmemo::HArray<IndexType> mDlg; //!< number of values in each column, size is mNumDiagonals
    hmemo::HArray<IndexType> mIlg; //!< number of values in each row, size is mNumRows
    hmemo::HArray<IndexType> mPerm; //!< position of each row in original matrix, size is mNumRows
    hmemo::HArray<IndexType> mJa; //!< column indices, size is mNumValues
    hmemo::HArray<ValueType> mValues; //!< non-zero values (+ optionally zeros in diagonal), size is mNumValues

private:

    /** Assigment with same storage type. */

    void assignJDS( const JDSStorage<ValueType>& other );

    /**
     *  Mandatory routine to check storage data for the diagonal property.
     */
    virtual bool checkDiagonalProperty() const;

    /**
     *  @brief Help routine for stable sort of ilg and setting up perm correctly
     *
     *  @param context is the preferred context where to execute it
     */
    void sortRows( hmemo::ContextPtr context );

    /** Help routine that sets up mDlg and allocates mJa and mValues after mIlg defined. */

    void setupData( hmemo::ContextPtr context );

    void print() const;

    SCAI_LOG_DECL_STATIC_LOGGER( logger ) //!< logger for this matrix format
};

} /* end namespace lama */

} /* end namespace scai */
