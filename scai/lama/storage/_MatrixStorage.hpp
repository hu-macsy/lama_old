/**
 * @file _MatrixStorage.hpp
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
 * @brief Definition of a common base class for all matrix storage classes.
 * @author Thomas Brandes
 * @date 27.04.2011
 */

#pragma once

// local library
#include <scai/lama/io/FileIO.hpp>
#include <scai/lama/storage/Format.hpp>

#include <scai/dmemo/Communicator.hpp>

// internal scai libraries
#include <scai/common/BinaryOp.hpp>
#include <scai/common/UnaryOp.hpp>
#include <scai/hmemo.hpp>

#include <scai/common/Factory.hpp>

// std
#include <ostream>
#include <vector>

namespace scai
{

// Forward declarations

namespace tasking
{
class SyncToken;
}

namespace dmemo
{
class Distribution;
class Halo;
class Redistributor;
}

namespace lama
{

/** Key type used for the _Matrix factory.
 *
 *  Note: own struct instead of std::pair to allow definition of operator <<
 */

typedef std::pair<Format, common::ScalarType> MatrixStorageCreateKeyType;

/** The class _MatrixStorage is the base class for all matrix storage classes
 supported by LAMA.

 In contratry to a matrix a matrix storage is never distributed. It contains
 local parts of a distributed matrix that are typically submatrices.

 This base class is also common for all supported value types for matrix
 elements.

 Note: Default copy constructor and assignment operator can be used.

 Each matrix storage has a context that will be taken as the preferred
 location where the matrix data resides and where operations on the matrix
 storage should take place.
 */

class COMMON_DLL_IMPORTEXPORT _MatrixStorage:
    public common::Factory<MatrixStorageCreateKeyType, _MatrixStorage*>,
    public common::Printable
{
public:

    /** Destructor. */

    virtual ~_MatrixStorage();

    /** This method returns the type name of a matrix storage, e.g. CSRStorage<double>.
     *
     *  This routine is more convenient for error messages than typeid(*this).name()
     */
    virtual const char* getTypeName() const = 0;

    /** Objects of this class should give info about their value type. */

    virtual common::ScalarType getValueType() const = 0;

    /** Clear the matrix storage, resets size to 0 x 0.
     *
     *  Clearing of matrix storage will not free allocated data, so resize to the original size
     *  is a cheap operation.
     */
    virtual void clear() = 0;

    /** Purge clears the matrix and even frees the allocated memory.
     */

    virtual void purge() = 0;

    /** Each storage provides a routine that prints the storage values. */

    virtual void print( std::ostream& stream ) const = 0;

    /** print vithout an argument for stream takes std::cout */

    void print() const
    {
        print( std::cout );
    }

    /** This method allocates new matrix storage for the matrix. The
     *  matrix contains only zero elements.
     *
     *  Note: This method contains an implicit clear, but not a purge.
     *        So allocation of a much smaller matrix compared to the existing
     *        one might result in some waste of memory.
     *
     *  Note: For dense matrices all values are set to zero; for sparse matrices
     *        zero values are usually not stored and no initialization is done.
     *
     *  Note: Operations are done on the context set for the matrix storage.
     */

    virtual void allocate( const IndexType numRows, const IndexType numColumns ) = 0;

    /** 
     *   This routine can be used to reset the number of columns. This might be necessary
     *   for IO where the number of columns is not explicitly specified and might be 
     *   adapted afterwards.
     *
     *   The default implementation assumes a sparse storage that is stored row-wise
     *   and so the size might be increased. 
     */
    virtual void resetNumColumns( const IndexType numColumns );

    /******************************************************************
     *  mContext  : setter / getter for preferred context              *
     ******************************************************************/

    /** @brief Set the preferred context for the matrix storage.
     *
     *  @param[in] context specifies where the storage should be allocated and where operations should be done
     */
    void setContextPtr( hmemo::ContextPtr context );

    /** @brief Getter for the preferred context of the storage data, returns pointer. */

    inline hmemo::ContextPtr getContextPtr() const;

    /** @brief Pure method that prefetches storage data into a given context.
     *
     *  @param context specifies location where data will resize
     */
    virtual void prefetch( const hmemo::ContextPtr context ) const = 0;

    /** @brief Method that prefetches storage data to its preferred location. */

    void prefetch() const
    {
        prefetch( mContext );
    }

    /** Will wait for all outstanding asynchronous data transfers. */

    virtual void wait() const = 0;

    /** Allow for additional row compression. */

    void setCompressThreshold( float ratio );

    /** Query the threshold  */

    float getCompressThreshold() const
    {
        return mCompressThreshold;
    }

    /** Getter routine for the number of rows. */

    inline IndexType getNumRows() const;

    /** Getter routine for the number of rows. */

    inline IndexType getNumColumns() const;

    virtual void writeAt( std::ostream& stream ) const;

    void resetDiagonalProperty();

    /** @brief Query if a matrix storage format has diagonal property, i.e. entries for each diagonal element.
     *
     *  Only if the diagonal property is given, the set/getDiagonal methods might be used to access it.
     *
     * The diagonal property for a dense matrix storage is always given but for a sparse matrix storage only
     * if the sparsity pattern includes all diagonal elements.
     *
     * Note for implementors: derived classes must set the flag mDiagonalProperty correctly.
     */
    inline bool hasDiagonalProperty() const;

    virtual Format getFormat() const = 0;

    /** This method sets storage for the identity matrix
     *
     *  @param[in] n is the size of the square matrix
     */

    virtual void setIdentity( const IndexType n ) = 0;

    /** This method resorts column indexes in such a way that the diagonal element is always the
     *  first one in a row.
     *
     *  This method throws an exception if the matrix storage is not square. Furthermore
     *  it throws an exception, if a diagonal element is zero, i.e. there is no entry for the diagonal
     *  element in a sparse format.
     */
    virtual void setDiagonalProperty();

    /** Get for each row the first column index with value entry.
     *  If diagonal flag is set, the column index will be the same as the (global) row
     *  index. I.e. for a local storage this routine gives the owned indexes to reconstruct
     *  the distribution.
     */

    virtual void getFirstColumnIndexes( hmemo::HArray<IndexType>& colIndexes ) const = 0;

    /******************************************************************
     *  General operations on a matrix                                 *
     ******************************************************************/

    /** This operation localizes the matrix storage of a full matrix to
     *  the part that is owned by this processor. This means that only
     *  the owned rows of the matrix will be kept.
     *
     *  Notes:
     *
     *    * routine can also be used if global is aliased with this matrix.
     *    * this routine is the same as an assign in case of a replicated distribution
     */

    virtual void localize( const _MatrixStorage& global, const dmemo::Distribution& rowDist );

    /** Get the total number of non-zero values in the matrix storage.
     *
     *  @return total number of non-zero values
     *
     *  An element is considered to be zero if its value is not stored explicitly
     *  (sparse format) or if its 'stored' value is explicitly zero. There is 
     *  the only exception that zero diagonal elements are also counted if 
     *  this->hasDiagonalProperty() is given.
     *
     *  This routine gives exactly the same number of elements that will be
     *  the size of the ja and values array when calling buildCSRData.
     *
     *  The default implementation uses the routine buildCSRSizes, but derived
     *  classes should override it if they can do it more efficiently.
     */
    virtual IndexType getNumValues() const;

    inline const hmemo::HArray<IndexType>& getRowIndexes() const;

    /**
     * @brief Get the number of entries in each row.
     *
     * @param[out] csrIA size array for rows, csrIA.size() == numRows
     *
     * If this->hasDiagonalProperty() is true, diagonal elements are also counted.
     *
     * \code
     *     DenseStorage<double> dense( ... )
     *     HArray<IndexType> ia;
     *     dense.buildCSRSizes( ia );
     * \endcode
     *
     * This routine must be provided by each matrix storage format.
     */

    virtual void buildCSRSizes( hmemo::HArray<IndexType>& csrIA ) const = 0;

    /**
     * @brief Get the matrix data of the storage in CSR format.
     *
     * @param[out] csrIA offset array for rows, csrIA.size() == numRows + 1
     * @param[out] csrJA column indexes, csrJA.size() == csrIA[ csrIA.size() ]
     * @param[out] csrValues are the non-zero matrix values, csrJA.size() == csrValues.size()
     *
     * The csr data will have the diagonal property if this->hasDiagonalProperty() is true.
     *
     * Note: This routine supports also type conversion between different value types.
     *
     * \code
     *     DenseStorage<double> dense( ... )
     *     HArray<IndexType> ia, ja;
     *     HArray<float> values;
     *     // get non-zero values of dense matrix as float values
     *     dense.buildCSRData( ia, ja, values );
     * \endcode
     *
     * This routine must be provided by each matrix storage format.
     */

    virtual void buildCSRData(
        hmemo::HArray<IndexType>& csrIA,
        hmemo::HArray<IndexType>& csrJA,
        hmemo::_HArray& csrValues ) const = 0;

    /** Each storage class must provide a routine to set CSR storage data.
     *
     *  @param numRows number of rows
     *  @param numColumns number of columns
     *  @param csrIA offset array ia for column indexes, size is numRows + 1
     *  @param csrJA are the column indexes of matrix entries, size is numValues
     *  @param csrValues are the values of matrix entries, size is numValues
     *
     *  \code
     *     _MatrixStorage& storage = ... 
     *     HArray<IndexType> ia( { 0, 1, 2, 4 } );
     *     HArray<IndexType> ja( { 0, 1, 2, 2 } );
     *     HArray<double> values( { 1.0, 1.0, 0.6, 0.4 } );
     *     storage.setCSRData( 3, 3, ia, ja, values );
     *  \endcode
     *
     *  This method checks, depending on the ASSERT_LEVEL, the input data whether
     *  it contains legal data.
     */
    virtual void setCSRData(
        const IndexType numRows,
        const IndexType numColumns,
        const hmemo::HArray<IndexType>& csrIA,
        const hmemo::HArray<IndexType>& csrJA,
        const hmemo::_HArray& csrValues ) = 0;

    /** Assign of matrix storage with any format or value type.
     *
     *  Format conversion and type conversion is done implicitly.
     *  Derived classes must implement this method.
     *
     *  @param[in] other is the matrix storage to be assigned
     *
     *  Note: a fallback solution is always to convert the other matrix to CSR format
     *        and assign it via buildCSRData and setCSRData.
     */
    virtual void assign( const _MatrixStorage& other ) = 0;

    /** The opposite routine to assign, for convenience as the other way around is
     *  sometimes more efficient
     */
    virtual void copyTo( _MatrixStorage& other ) const = 0;

    /**
     *  Extract a contiguous block of n rows in a new storage
     *
     *  @param[out] other  will contain the corresponding storage of n x mNumColumns
     *  @param[in]  first  index of first row to extract
     *  @param[in]  n      number of rows to extract
     */
    virtual void copyBlockTo( _MatrixStorage& other, const IndexType first, const IndexType n ) const = 0;

    /** Override default assignment operator. */

    _MatrixStorage& operator=( const _MatrixStorage& other );

    /** Delete the move assignment operator in this base class as there is no virtual move semantic
     *  for the derived classes. 
     */

    _MatrixStorage& operator=( _MatrixStorage&& other ) = delete;

    /******************************************************************
     *   Help routines (ToDo: -> HArrayUtils ?? )                   *
     ******************************************************************/

    /** Help routines to convert arrays with sizes to offsets and vice versa */

    static void offsets2sizes( hmemo::HArray<IndexType>& offsets );

    static void offsets2sizes( hmemo::HArray<IndexType>& sizes, const hmemo::HArray<IndexType>& offsets );

    static IndexType sizes2offsets( hmemo::HArray<IndexType>& offsets );

    static IndexType sizes2offsets( hmemo::HArray<IndexType>& offsets, const hmemo::HArray<IndexType>& sizes, const hmemo::ContextPtr loc );

    /** Returns the number of bytes needed for the current matrix.
     *
     *  Note: This routine does not tell how many memory is really allocated. Storage data
     *        might be allocated on more than one device. Furthermore, it is possible that
     *        arrays have more memory reserved than needed for its current size.
     */

    /******************************************************************
     *   Query routines                                                *
     ******************************************************************/

    size_t getMemoryUsage() const;

    virtual void check( const char* msg ) const = 0;

    /** Each matrix storage must provide a routine that makes a new copy
     *  of the input matrix (same format and same value type).
     *
     *  The implementations in derived classes should use covariant
     *  return types to allow polymorphic creation of objects.
     *
     *  Note: as shared_ptr cannot be used for for covariant return types
     *        this routine returns directly a new created object that should
     *        be wrapped as a shared pointer at calling site.
     */

    virtual _MatrixStorage* copy() const = 0;

    /**  Each matrix storage must provide a routine that creates a new zero matrix storage
      *  with same attributes (size, type, context).
      *
      *  \code
      *      const _MatrixStorage& anyA = ...
      *      std::unique_ptr<_MatrixStorage> tmpA( anyA.newMatrixStorage() );
      *  \endcode
      *
      *  Note: this method returns a pointer to a new allocated object and it is recommended
      *        to keep it as unique or shared pointer to guarantee that it is freed.
      */
    virtual _MatrixStorage* newMatrixStorage( const IndexType numRows, const IndexType numColumns ) const = 0;

    /**
     *  Calling newMatrixStorage without arguments returns storage with same size as this storage.
     */
    _MatrixStorage* newMatrixStorage() const
    {
        return newMatrixStorage( getNumRows(), getNumColumns() );
    }

    /**  Returning a pair of storage format and value type for using with
      *  the factory
      */

    virtual MatrixStorageCreateKeyType getCreateValue() const = 0;

    /**
     * @brief read the matrix storage from an input file
     *
     * @param[in] fileName is the name of the input file (suffix must be added according to the file type)
     * @param[in] firstRow is the first row to read
     * @param[in] nRows    specifies the number of rows to read, defaults to number of rows of full storage - firstRow
     *
     * Note: default argument for nRows is invalidIndex as the number of rows in full storage might not be known
     */
    virtual void readFromFile(
        const std::string& fileName,
        const IndexType firstRow = 0,
        const IndexType nRows = invalidIndex ) = 0;

    /**
     * @brief write the matrix storage to an output file
     *
     * @param[in] fileName  is the name of the output file (suffix must be added according to the file type)
     * @param[in] type      format of the output file ("frm" for SAMG, "mtx" for _MatrixMarket), default is to decide by suffix
     * @param[in] dataType  representation type for output values, default is same type as matrix values
     * @param[in] indexType representation type for row/col index values (default is settings of FileIO)
     * @param[in] fileMode  use BINARY or FORMATTED to force a certain mode, otherwise DEFAULT
     *
     * If one of the arguments dataType, indexType or fileMode is set, it will overwrite
     * any setting specified by the corresponding environment variables SCAI_IO_TYPE_DATA, SCAI_IO_TYPE_INDEX
     * or SCAI_IO_BINARY
     */

    virtual void writeToFile(
        const std::string& fileName,
        const std::string& type = "",
        const common::ScalarType dataType = common::ScalarType::UNKNOWN,
        const common::ScalarType indexType = common::ScalarType::UNKNOWN,
        const FileIO::FileMode fileMode = FileIO::DEFAULT_MODE  ) const = 0;

    virtual bool checkSymmetry() const = 0;

protected:

    /** Constructor with specification of context
     * 
     *  @param[in] numRows     number of rows the matrix storage has
     *  @param[in] numColumns  number of columns the matrix storage has
     *  @param[in] ctx         will be the context of the storage
     *
     *  The context decides where memory will be allocated and operations on the storage are executed.
     */
    _MatrixStorage( IndexType numRows, IndexType numColumns, hmemo::ContextPtr ctx  );

    /** Move constructor that might be used by derived classes */

    _MatrixStorage( _MatrixStorage&& other ) noexcept;

    _MatrixStorage( const _MatrixStorage& other );

    /** Swaps this with other.
     *
     * swap is protected to avoid accidently wrong swaps of base classes which do not
     * implement their own swap.
     *
     * @param[in,out] other the _MatrixStorage to swap this with
     */
    void swap( _MatrixStorage& other );

    /** Set the dimension of the storage.
     *
     *  @param numRows, numColumns specifiy the size of the matrix
     */
    void setDimension( const IndexType numRows, const IndexType numColumns );

    /** 
     *  @brief Splitting up a matrix storage returns sizes and sets this storage to zero storage.
     *
     *  @param[out] numRows same as getNumRows()
     *  @param[out] numColumns same as getNumColumns()
     *
     *  This storage itself becomes a zero storage.
     */
    void splitUp( IndexType& numRows, IndexType& numColumns );

    void _assignTranspose( const _MatrixStorage& other );

    /** This methods provides the assignment operator semantic for the member variables of this class only. */

    void _assign( const _MatrixStorage& other );

    /** This method provides the move semantic for member variables of this base class only. */

    void moveImpl( _MatrixStorage&& other );

    /** Returns the number of bytes needed for the current matrix.
     *
     *  This pure function must be implemented by each derived class. Relevant for the number
     *  of bytes is the current size of the used (LAMA) arrays.
     *
     *  Note: This routine does not tell how many memory is really allocated. Storage data
     *        might be allocated on more than one device. Furthermore, it is possible that
     *        arrays have more memory reserved than needed for its current size.
     */

    virtual size_t getMemoryUsageImpl() const = 0;

private:

    IndexType mNumRows; //!< numbers of matrix rows

    IndexType mNumColumns; //!< number of matrix columns

protected:

    hmemo::HArray<IndexType> mRowIndexes; //!< used in case of sparse representation of ia

    float mCompressThreshold; //!< ratio at which compression is done, 0 for never, 1 for always

    bool mDiagonalProperty; //!< if true, diagonal elements are always stored at first position in each row

    hmemo::ContextPtr mContext;//!< preferred context for the storage

    SCAI_LOG_DECL_STATIC_LOGGER( logger ); //!< logger for this matrix format

protected:

    /** checkDiagonalProperty checks if the diagonal property of this is full filled.
     *
     * @return true if the diagonal property is fulfilled for the matrix data
     */

    virtual bool checkDiagonalProperty() const = 0;

};

/* ------------------------------------------------------------------------- */
/* Implementation of inline methods                                          */
/* ------------------------------------------------------------------------- */

hmemo::ContextPtr _MatrixStorage::getContextPtr() const
{
    return mContext;
}

IndexType _MatrixStorage::getNumRows() const
{
    return mNumRows;
}

IndexType _MatrixStorage::getNumColumns() const
{
    return mNumColumns;
}

bool _MatrixStorage::hasDiagonalProperty() const
{
    return mDiagonalProperty;
}

const hmemo::HArray<IndexType>& _MatrixStorage::getRowIndexes() const
{
    return mRowIndexes;
}

} /* end namespace lama */

} /* end namespace scai */
