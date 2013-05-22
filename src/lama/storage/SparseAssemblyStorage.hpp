/**
 * @file SparseAssemblyStorage.hpp
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
 * @brief SparseAssemblyStorage.hpp
 * @author Jiri Kraus
 * @date 07.11.2011
 * @since 1.0.0
 */
#ifndef LAMA_SPARSEASSEMBLYSTORAGE_HPP_
#define LAMA_SPARSEASSEMBLYSTORAGE_HPP_

// for dll_import
#include <lama/config.hpp>

// others
#include <lama/storage/CRTPMatrixStorage.hpp>

namespace lama
{

/**
 * @brief The template SparseAssemblyStorage is a convenient intermediate container to construct
 *        a sparse matrix which structure is not known a priori.
 *
 * The SparseAssemblyStorage should allow a convenient construction of a sparse matrix if the
 * structure of this in not know prior to its construction. The complexity of inserting a new
 * element is constant (amortized time, reallocation may happen).
 *
 * Due to the use of std::vector, this storage will only reside on the host.
 *
 *  @tparam T is the value type of the matrix values.
 */
template<typename T>
class LAMA_DLL_IMPORTEXPORT SparseAssemblyStorage:

    public CRTPMatrixStorage<SparseAssemblyStorage<T>,T>
{
public:
    typedef T ValueType; //!< This is the type of the matrix values.

    /** Getter of the type name of the matrix storage format. */

    static const char* typeName();

    /**
     * @brief The default constructor creates an empty SparseAssemblyStorage.
     */
    SparseAssemblyStorage();

    /**
     * @brief The copy constructor creates a SparseAssemblyStorage which is a
     *        copy of the passed SparseAssemblyStorage.
     *
     * @param[in]   other the SparseAssemblyStorage to take a copy from.
     */
    SparseAssemblyStorage( const SparseAssemblyStorage<ValueType>& other );

    /**
     * @brief This constructor creates a SparseAssemblyStorage which is
     *        initialized with the values of the passed MatrixStorage.
     *
     * @param[in]   other the MatrixStorage to take the values for initialization.
     */
    explicit SparseAssemblyStorage( const _MatrixStorage& other );

    /**
     * @brief This constructor creates a SparseAssemblyStorage with the given
     *        number of rows and columns.
     *
     * This constructor creates a SparseAssemblyStorage with the given number of
     * rows and columns. As a performance optimization the user can specify the
     * optional parameter numValuesPerRow, which determines the approximate number of
     * none zeros per row needed.
     *
     * @param[in] numRows           the number of rows the matrix should have.
     * @param[in] numColumns        the number of columns the matrix should have.
     * @param[in] numValuesPerRow   the storage that should be reserver for each row.
     *                              (Default: 10)
     */
    SparseAssemblyStorage( const IndexType numRows, const IndexType numColumns, const IndexType numValuesPerRow = 10 );

    /**
     * @brief The destructor releases the allocated memory.
     */
    virtual ~SparseAssemblyStorage();

    virtual void writeAt( std::ostream& stream ) const;

    /**
     * @brief Implementation of pure method for _MatrixStorage.
     */

    virtual void clear();

    /** Implementation of pure method of class MatrixStorage. */

    virtual void purge();

    /** Swaps this storage with other storage of the same type.
     *
     *  @param[in,out] other is the storage to be swapped with this one
     */
    void swap( SparseAssemblyStorage<ValueType>& other );

    /** Implementation of pure method _MatrixStorage::allocate */

    void allocate( const IndexType numRows, const IndexType numColumns );

    /**
     * @brief The assignment operator sets the values of this to the values of
     *        the passed SparseAssemblyStorage
     *
     * @param[in]   other the SparseAssemblyStorage to take a copy from.
     * @return      a reference to this.
     */
    SparseAssemblyStorage<T>& operator=( const SparseAssemblyStorage<T>& other );

    /** Assignment operator for arbitrary storages. */

    SparseAssemblyStorage<T>& operator=( const _MatrixStorage& other );

    /**
     * @brief the access operator returns a copy of the element at the
     *        queried position.
     *
     * The access operator returns a copy of the element at the queried
     * position. The runtime is linear with respect to the number of none zero
     * values in the queried row, because the rows do not need to be sorted.
     * A zero is returned if the there is no element at the queried position.
     *
     * @param[in] i the row to get the value at.
     * @param[in] j the column to get the value at.
     * @return      a copy of the value at (i,j) or zero.
     */
    ValueType operator()( const IndexType i, const IndexType j ) const;

    /** Implementation for MatrixStorage::maxNorm */

    virtual ValueType maxNorm() const;

    ValueType getValue( const IndexType i, const IndexType j ) const
    {
        return ( *this )( i, j );
    }

    /** Getter routine for the enum value that stands for this format. */

    virtual MatrixStorageFormat getFormat() const;

    /** Implementation of pure method.  */

    virtual void setIdentity( const IndexType size );

    /**
     * @brief sets the element at the queried position.
     *
     * sets the element at the queried position. The runtime is linear with
     * respect to the number of none zero values in the queried row, because the
     * rows do not need to be sorted. Although inserting an element
     * at the end has amortized constant complexity a reallocation may be needed
     * if there is not enough free storage is available in the queried row.
     *
     * @param[in] i         the row to get the value at.
     * @param[in] j         the column to get the value at.
     * @param[in] value     the value of the at (i,j)
     */
    void set( const IndexType i, const IndexType j, const ValueType value );

    /**
     * @brief getJa returns a constant reference to the column index array of the
     *        given row.
     *
     * @param[in] i the row to get the column index array for.
     * @return      a const reference to the column index array of row i.
     */
    const std::vector<IndexType>& getJa( const IndexType i ) const;

    std::vector<IndexType>& getJa( const IndexType i );

    /**
     * @brief getValues returns a constant reference to the values array of the
     *        given row.
     *
     * @param[in] i the row to get the value array for.
     * @return      a const reference to the value array of row i.
     */
    const std::vector<ValueType>& getValues( const IndexType i ) const;

    /** Initiate an asynchronous data transfer to a specified location. */

    virtual void prefetch( const ContextPtr location ) const;

    /** Will wait for all outstanding asynchronous data transfers. */

    void wait() const;

    /**
     * @brief getNumValues returns the number values stored in this matrix.
     *
     * @return  the number of values stored in this matrix.
     */
    IndexType getNumValues() const;

    /** Implementation for pure method is provided. */

    virtual size_t getMemoryUsageImpl() const;

    /**
     * @brief setRow sets the given row to the passed values.
     *
     * setRow sets the given row to the passed values. The values are copied
     * into a possibly reallocated internal row.
     *
     * @param[in] i     the row to set.
     * @param[in] ja    the column vector of row i.
     * @param[in] a     the values of row i.
     */
    void setRow( const IndexType i, const LAMAArrayConstView<IndexType>& ja, const LAMAArrayConstView<ValueType>& a );

    /**
     * @brief fixDiagonalProperty fixes the diagonal property of row.
     *
     * fixDiagonalProperty fixes the diagonal property of row. If needed
     * a zero entry is added.
     *
     * @param[in] i     the row to fix the diagonal property.
     */
    void fixDiagonalProperty( const IndexType i );

    /**
     * @brief setNumColumns sets the number of columns.
     *
     * setNumColumns sets the number of columns. This should be done only during
     * the assembly of a matrix, not after diagonal properties have been fixed!
     *
     * @param[in] i     the number of columns.
     */
    void setNumColumns( const IndexType i );

    /** Template method for getting row. */

    template<typename OtherType>
    void getRowImpl( LAMAArray<OtherType>& row, const IndexType i ) const;

    /** Typed version of getDiagonal
     *
     * @param[out] diagonal is the typed LAMA array for output
     */
    template<typename OtherValueType>
    void getDiagonalImpl( LAMAArray<OtherValueType>& diagonal ) const;

    /** Typed version of setDiagonal
     *
     * @param[in] diagonal is the typed LAMA array for input
     */
    template<typename OtherValueType>
    void setDiagonalImpl( const LAMAArray<OtherValueType>& diagonal );

    /** Implementation of pure method. */

    void setDiagonalImpl( const Scalar scalar );

    /******************************************************************
     *  Scaling of elements in a matrix                                *
     ******************************************************************/

    /** Template version used for virtual routine scale with known value type. */

    template<typename OtherType>
    void scaleImpl( const LAMAArray<OtherType>& values );

    /** Implementation of pure method.  */

    void scaleImpl( const Scalar value );

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
        const ContextPtr loc ) const;

    /**
     * @brief fills Assembly sparse matrix by csr sparse data.
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

    /** Test the storage data for inconsistencies.
     *
     *  @param[in] msg is additional message string that should be used to identify calling routine
     *  @throws Exception in case of any inconsistency.
     */
    void check( const char* msg ) const;

    virtual SparseAssemblyStorage* copy() const
    {
        return new SparseAssemblyStorage( *this );
    }

    /** Implementation of MatrixStorage::create for derived class. */

    virtual SparseAssemblyStorage* create() const
    {
        return new SparseAssemblyStorage();
    }

    using MatrixStorage<ValueType>::assign;

private:

    void print() const;

    using CRTPMatrixStorage<SparseAssemblyStorage<T>,T>::mDiagonalProperty;
    using CRTPMatrixStorage<SparseAssemblyStorage<T>,T>::mNumRows;
    using CRTPMatrixStorage<SparseAssemblyStorage<T>,T>::mNumColumns;

    struct Row
    {
        std::vector<IndexType> ja;
        std::vector<ValueType> values;
        Row();
        Row( const typename SparseAssemblyStorage<T>::Row& other );
        Row( const IndexType numValuesPerRow );
        Row& operator=( const Row& other );
        void reserve( const IndexType numValuesPerRow );
        void scale( const ValueType val );
    };

    std::vector<Row> mRows;

    IndexType mNumValues; //!< count number of all row elements

    /**
     * @brief checks if in each row the diagonal element is stored first.
     *
     */
    virtual bool checkDiagonalProperty() const;

    LAMA_LOG_DECL_STATIC_LOGGER( logger )
};

} //namespace lama

#endif // LAMA_SPARSEASSEMBLYSTORAGE_HPP_
