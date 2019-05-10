/**
 * @file AssemblyStorage.hpp
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
 * @brief AssemblyStorage.hpp
 * @author Thomas, Jiri Kraus
 * @date 07.05.2019
 */
#pragma once

// for dll_import
#include <scai/common/config.hpp>

// base classes
#include <scai/lama/storage/MatrixStorage.hpp>
#include <scai/lama/mepr/StorageWrapper.hpp>

namespace scai
{

namespace lama
{

/**
 * @brief The template AssemblyStorage is a convenient intermediate container to construct
 *        a sparse matrix which structure is not known a priori.
 *
 * The AssemblyStorage should allow a convenient construction of a sparse matrix if the
 * structure of this in not know prior to its construction. The complexity of inserting a new
 * element is constant (amortized time, reallocation may happen).
 *
 * Due to the use of std::vector, this storage will only reside on the host.
 *
 *  @tparam ValueType is the value type of the matrix values.
 */
template<typename ValueType>
class COMMON_DLL_IMPORTEXPORT AssemblyStorage:

    public MatrixStorage<ValueType>,
    public _MatrixStorage::Register<AssemblyStorage<ValueType> >    // register at factory

{
public:

    /** Getter of the type name of the matrix storage format. */

    static const char* typeName();

    /** Implementation of pure method _MatrixStorage:getTypeName    */

    virtual const char* getTypeName() const;

    /** Implementation of pure method _MatrixStorage:getCreateValue    */

    virtual MatrixStorageCreateKeyType getCreateValue() const;

    /**
     * @brief The default constructor creates an empty AssemblyStorage.
     */
    AssemblyStorage();

    /**
     * @brief The copy constructor creates a AssemblyStorage which is a
     *        copy of the passed AssemblyStorage.
     *
     * @param[in]   other the AssemblyStorage to take a copy from.
     */
    AssemblyStorage( const AssemblyStorage<ValueType>& other );

    /**
     * @brief This constructor creates a AssemblyStorage which is
     *        initialized with the values of the passed MatrixStorage.
     *
     * @param[in]   other the MatrixStorage to take the values for initialization.
     */
    explicit AssemblyStorage( const _MatrixStorage& other );

    /**
     * @brief This constructor creates a AssemAssembith the given
     *        number of rows and columns.
     *
     * This constructor creates a AssemblyStorage with the given number of
     * rows and columns. As a performance optimization the user can specify the
     * optional parameter numValuesPerRow, which determines the approximate number of
     * none zeros per row needed.
     *
     * @param[in] numRows           the number of rows the matrix should have.
     * @param[in] numColumns        the number of columns the matrix should have.
     * @param[in] numValuesPerRow   the storage that should be reserver for each row.
     *                              (Default: 10)
     */
    AssemblyStorage( const IndexType numRows, const IndexType numColumns, const IndexType numValuesPerRow = 10 );

    /**
     * @brief The destructor releases the allocated memory.
     */
    virtual ~AssemblyStorage();

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
    void swapImpl( AssemblyStorage<ValueType>& other );

    /** Implementation for _MatrixStorage::swap */

    virtual void swap( _MatrixStorage& other );

    /** Implementation of pure method _MatrixStorage::allocate */

    void allocate( const IndexType numRows, const IndexType numColumns );

    /**
     * @brief The assignment operator sets the values of this to the values of
     *        the passed AssemblyStorage
     *
     * @param[in]   other the AssemblyStorage to take a copy from.
     * @return      a reference to this.
     */
    AssemblyStorage<ValueType>& operator=( const AssemblyStorage<ValueType>& other );

    /** Assignment operator for arbitrary storages. */

    AssemblyStorage<ValueType>& operator=( const _MatrixStorage& other );

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

    virtual RealType<ValueType> l1Norm() const;

    /** Implementation for MatrixStorage::l2Norm */

    virtual RealType<ValueType> l2Norm() const;

    /** Implementation for MatrixStorage::maxNorm */

    virtual RealType<ValueType> maxNorm() const;

    ValueType getValue( const IndexType i, const IndexType j ) const
    {
        return ( *this )( i, j );
    }

    /** Implementation of pure method MatrixStorage<ValueType>::setValue for ELL storage */

    void setValue( const IndexType i, const IndexType j, const ValueType val,
                   const common::BinaryOp = common::BinaryOp::COPY );

    /** Getter routine for the enum value that stands for this format. */

    virtual Format getFormat() const;

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

    /** Implementation of pure method MatrixStorage<ValueType>::getRow */

    virtual void getRow( hmemo::HArray<ValueType>& row, const IndexType i ) const;

    /** Implementation of pure method MatrixStorage<ValueType>::getSparseRow */

    virtual void getSparseRow( hmemo::HArray<IndexType>& jA, hmemo::HArray<ValueType>& values, const IndexType i ) const;

    /** Implementation of pure method MatrixStorage<ValueType>::setRow */

    virtual void setRow( const hmemo::HArray<ValueType>& row, const IndexType i, const common::BinaryOp op );

    /** Implementation of pure method MatrixStorage<ValueType>::getColumn */

    virtual void getColumn( hmemo::HArray<ValueType>& column, const IndexType j ) const;

    /** Implementation of pure method MatrixStorage::getSparseColumn */

    virtual void getSparseColumn( hmemo::HArray<IndexType>& iA, hmemo::HArray<ValueType>& values, const IndexType j ) const;

    /** Implementation of pure method MatrixStorage<ValueType>::setColumn */

    virtual void setColumn( const hmemo::HArray<ValueType>& column, const IndexType j, const common::BinaryOp op );

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

    /** Template version used for virtual routine scale with known value type. */

    template<typename OtherType>
    void scaleImpl( const hmemo::HArray<OtherType>& values );

    /** Implementation of pure method.  */

    void scaleImpl( const ValueType value );

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

    /** Implementation for _MatrixStorage::buildCSRSizes */

    void buildCSRSizes( hmemo::HArray<IndexType>& ia ) const;

    /** Implementation for _MatrixStorage::buildCSRData */

    void buildCSRData( hmemo::HArray<IndexType>& csrIA, hmemo::HArray<IndexType>& csrJA, hmemo::_HArray& csrValues ) const;

    /**
     * @brief fills assembly matrix by csr sparse data.
     *
     * @param[in] numRows    number of rows
     * @param[in] numColumns the number of columsns
     * @param[in] ia         row offset array for access to ja and values
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

    /** Implementation of pure method _MatrixStorage::setCSRData for this class.  */

    void setCSRData(
        const IndexType numRows,
        const IndexType numColumns,
        const hmemo::HArray<IndexType>& ia,
        const hmemo::HArray<IndexType>& ja,
        const hmemo::_HArray& values );

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

    virtual AssemblyStorage* newMatrixStorage( const IndexType numRows, const IndexType numColumns ) const;

    /** Implementation of MatrixStorage::copy for derived class. */
    virtual AssemblyStorage* copy() const
    {
        return new AssemblyStorage( *this );
    }

    /**
     * @brief Implemenation of pure method MatrixStorage<ValueType>::assignDiagonal
     */
    virtual void assignDiagonal( const hmemo::HArray<ValueType>& diagonal );

    /**
     * @brief Implementation of pure method _MatrixStorage::assign
     */
    virtual void assign( const _MatrixStorage& other );

    /**
     *  @brief Implemenation of assignments for this class
     */
    template<typename OtherValueType>
    void assignImpl( const MatrixStorage<OtherValueType>& other );

    using _MatrixStorage::prefetch;
    using _MatrixStorage::getContextPtr;

    using _MatrixStorage::getNumRows;
    using _MatrixStorage::getNumColumns;

    using MatrixStorage<ValueType>::getValueType;

    /** Implementation of MatrixStorage::matrixTimesVector for Assembly */

    virtual void matrixTimesVector(
        hmemo::HArray<ValueType>& result,
        const ValueType alpha,
        const hmemo::HArray<ValueType>& x,
        const ValueType beta,
        const hmemo::HArray<ValueType>& y,
        const common::MatrixOp op ) const;

private:

    /** Implementation of pure method. */

    void print( std::ostream& stream ) const;

    struct Row
    {
        std::vector<IndexType> ja;
        std::vector<ValueType> values;
        Row();
        Row( const typename AssemblyStorage<ValueType>::Row& other );
        Row( const IndexType numValuesPerRow );
        Row& operator=( const Row& other );
        void reserve( const IndexType numValuesPerRow );
        void scale( const ValueType val );
        void conj();
    };

    std::vector<Row> mRows;

    IndexType mNumValues; //!< count number of all row elements

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
