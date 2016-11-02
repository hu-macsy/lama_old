/**
 * @file SparseAssemblyStorage.hpp
 *
 * @license
 * Copyright (c) 2009-2016
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
 * @brief SparseAssemblyStorage.hpp
 * @author Jiri Kraus
 * @date 07.11.2011
 */
#pragma once

// for dll_import
#include <scai/common/config.hpp>

// base classes
#include <scai/lama/storage/CRTPMatrixStorage.hpp>

namespace scai
{

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
 *  @tparam ValueType is the value type of the matrix values.
 */
template<typename ValueType>
class COMMON_DLL_IMPORTEXPORT SparseAssemblyStorage:

    public CRTPMatrixStorage<SparseAssemblyStorage<ValueType>, ValueType>,
    public _MatrixStorage::Register<SparseAssemblyStorage<ValueType> >    // register at factory
{
public:

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
    void swapImpl( SparseAssemblyStorage<ValueType>& other );

    /** Implementation for _MatrixStorage::swap */

    virtual void swap( _MatrixStorage& other );

    /** Implementation of pure method _MatrixStorage::allocate */

    void allocate( const IndexType numRows, const IndexType numColumns );

    /**
     * @brief The assignment operator sets the values of this to the values of
     *        the passed SparseAssemblyStorage
     *
     * @param[in]   other the SparseAssemblyStorage to take a copy from.
     * @return      a reference to this.
     */
    SparseAssemblyStorage<ValueType>& operator=( const SparseAssemblyStorage<ValueType>& other );

    /** Assignment operator for arbitrary storages. */

    SparseAssemblyStorage<ValueType>& operator=( const _MatrixStorage& other );

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

    /** Implementation for MatrixStorage::l1Norm */

    virtual ValueType l1Norm() const;

    /** Implementation for MatrixStorage::l2Norm */

    virtual ValueType l2Norm() const;

    /** Implementation for MatrixStorage::maxNorm */

    virtual ValueType maxNorm() const;

    ValueType getValue( const IndexType i, const IndexType j ) const
    {
        return ( *this )( i, j );
    }

    /** Implementation of pure method MatrixStorage<ValueType>::setValue for ELL storage */

    void setValue( const IndexType i, const IndexType j, const ValueType val,
                   const utilskernel::binary::BinaryOp op = utilskernel::binary::COPY );

    /** Getter routine for the enum value that stands for this format. */

    virtual Format::MatrixStorageFormat getFormat() const;

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

    virtual void prefetch( const hmemo::ContextPtr location ) const;

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
    void setSparseRow( const IndexType i, const hmemo::HArray<IndexType>& ja, const hmemo::HArray<ValueType>& a );

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

    /** Template version of getRow */

    template<typename OtherType>
    void getRowImpl( hmemo::HArray<OtherType>& row, const IndexType i ) const;

    /** Template version of setRow */

    template<typename OtherType>
    void setRowImpl( const hmemo::HArray<OtherType>& row, const IndexType i,
                     const utilskernel::binary::BinaryOp op );

    /** Template version of getColumn */

    template<typename OtherType>
    void getColumnImpl( hmemo::HArray<OtherType>& column, const IndexType j ) const;

    /** Template version of setColumn */

    template<typename OtherType>
    void setColumnImpl( const hmemo::HArray<OtherType>& column, const IndexType j,
                        const utilskernel::binary::BinaryOp op );

    /** Typed version of getDiagonal
     *
     * @param[out] diagonal is the typed LAMA array for output
     */
    template<typename OtherValueType>
    void getDiagonalImpl( hmemo::HArray<OtherValueType>& diagonal ) const;

    /** Typed version of setDiagonal
     *
     * @param[in] diagonal is the typed LAMA array for input
     */
    template<typename OtherValueType>
    void setDiagonalImpl( const hmemo::HArray<OtherValueType>& diagonal );

    /** Implementation of pure method. */

    void setDiagonalImpl( const ValueType value );

    /******************************************************************
     *  Scaling of elements in a matrix                                *
     ******************************************************************/

    /** Template version used for virtual routine scale with known value type. */

    template<typename OtherType>
    void scaleImpl( const hmemo::HArray<OtherType>& values );

    /** Implementation of pure method.  */

    void scaleImpl( const ValueType value );

    /** Implementation of pure method.  */

    void conj();

    /** General routine to build any kind of CSR storage.
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
        const hmemo::HArray<IndexType>& ia,
        const hmemo::HArray<IndexType>& ja,
        const hmemo::HArray<OtherValueType>& values,
        const hmemo::ContextPtr loc );

    /**
     * @brief fills Assembly sparse matrix by dia sparse data.
     *
     * @param[in] numRows      number of rows
     * @param[in] numColumns   number of columns
     * @param[in] numDiagonals the number of stored diagonals
     * @param[in] offsets      raw pointer of the input csr sparse matrix
     * @param[in] values       the data values of the input csr sparse matrix
     * @param[in] loc          is the context where filling takes place
     */
    template<typename OtherValueType>
    void setDIADataImpl(
        const IndexType numRows,
        const IndexType numColumns,
        const IndexType numDiagonals,
        const hmemo::HArray<IndexType>& offsets,
        const hmemo::HArray<OtherValueType>& values,
        const hmemo::ContextPtr loc ) __attribute__( ( noinline ) );

    /** Test the storage data for inconsistencies.
     *
     *  @param[in] msg is additional message string that should be used to identify calling routine
     *  @throws Exception in case of any inconsistency.
     */
    void check( const char* msg ) const;

    /** Implementation of MatrixStorage::newMatrixStorage for derived class. */
    virtual SparseAssemblyStorage* newMatrixStorage() const;

    /** Implementation of MatrixStorage::copy for derived class. */
    virtual SparseAssemblyStorage* copy() const
    {
        return new SparseAssemblyStorage( *this );
    }

    using MatrixStorage<ValueType>::assign;

private:

    /** Implementation of pure method. */

    void print( std::ostream& stream ) const;

    using CRTPMatrixStorage<SparseAssemblyStorage<ValueType>, ValueType>::mDiagonalProperty;
    using CRTPMatrixStorage<SparseAssemblyStorage<ValueType>, ValueType>::mNumRows;
    using CRTPMatrixStorage<SparseAssemblyStorage<ValueType>, ValueType>::mNumColumns;

    struct Row
    {
        std::vector<IndexType> ja;
        std::vector<ValueType> values;
        Row();
        Row( const typename SparseAssemblyStorage<ValueType>::Row& other );
        Row( const IndexType numValuesPerRow );
        Row& operator=( const Row& other );
        void reserve( const IndexType numValuesPerRow );
        void scale( const ValueType val );
        void conj();
    };

    std::vector<Row> mRows;

    IndexType mNumValues; //!< count number of all row elements

    /**
     * @brief checks if in each row the diagonal element is stored first.
     *
     */
    virtual bool checkDiagonalProperty() const;

    SCAI_LOG_DECL_STATIC_LOGGER( logger )

    static std::string initTypeName();

public:

    // static create method that will be used to register at MatrixStorage factory

    static _MatrixStorage* create();

    // key for factory

    static MatrixStorageCreateKeyType createValue();
};

} /* end namespace lama */

} /* end namespace scai */
