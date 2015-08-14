/**
 * @file MatrixStorage.hpp
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
 * @brief Definition of a common base class for all matrix storage formats.
 * @author Thomas Brandes
 * @date 27.04.2011
 * @since 1.0.0
 */

#pragma once

#include <ostream>

// others
#include <scai/memory.hpp>
#include <scai/lama/Communicator.hpp>
#include <scai/lama/Scalar.hpp>

#include <scai/lama/io/FileType.hpp>

namespace tasking
{
class SyncToken;
}

using memory::Context;
using memory::ContextPtr;
using memory::ContextArray;
using memory::LAMAArray;
using memory::LAMAArrayRef;

namespace scai
{

namespace lama
{

// Forward declaration

class Distribution;

template<typename ValueType> class CSRStorage;

/** Enumeration type for different matrix storage formats.
 *
 *  Note: operator<< is implemented for this type and should be adapted in case of modifications.
 */
namespace Format
{
typedef enum
{
    CSR, //!< Compressed Sparse Row
    ELL, //!< ELLPack
    DIA, //!< Diagonal
    JDS, //!< Jagged Diagonal Storage
    COO, //!< Coordinate list
    DENSE, //!< Dense, all elements are stored
    ASSEMBLY, //!<  Matrix storage used for assembling of values
    UNDEFINED //!<  Default value
} MatrixStorageFormat;
}
;

using Format::MatrixStorageFormat;
//!< useful abbreviation

COMMON_DLL_IMPORTEXPORT const char* format2Str( const MatrixStorageFormat storageFormat );

COMMON_DLL_IMPORTEXPORT MatrixStorageFormat str2Format( const char* str );

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

class COMMON_DLL_IMPORTEXPORT _MatrixStorage: public Printable
{
public:

    /** Default constructor */

    _MatrixStorage();

    /** Destructor. */

    virtual ~_MatrixStorage();

    /** This method returns the type name of a matrix storage, e.g. CSRStorage<double>.
     *
     *  This routine is more convenient for error messages than typeid(*this).name()
     */

    virtual const char* getTypeName() const = 0;

    /** Set the dimension of the storage.
     *
     *  @param numRows, numColumns specifiy the size of the matrix
     */

    void setDimension( const IndexType numRows, const IndexType numColumns );

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

    /** This method allocates new matrix storage for the matrix. The
     *  matrix contains only zero elements.
     *
     *  Note: This method contains an implicit clear, but not a purge.
     *        So allocation of a much smaller matrix compared to the existing
     *        one might result in some waste of memory.
     */

    virtual void allocate( const IndexType numRows, const IndexType numColumns ) = 0;

    /******************************************************************
     *  mContext  : setter / getter for preferred context              *
     ******************************************************************/

    /** @brief Set the preferred context for the matrix storage.
     *
     *  @param[in] context specifies where the storage should be allocated and where operations should be done
     */

    void setContext( ContextPtr context );

    /** @brief Getter for the preferred context of the storage data, returns pointer. */

    inline ContextPtr getContextPtr() const;

    /** @brief Getter for the preferred context of the storage data. */

    inline const Context& getContext() const;

    /** @brief Pure method that prefetches storage data into a given context.
     *
     *  @param context specifies location where data will resize
     */

    virtual void prefetch( const ContextPtr context ) const = 0;

    /** @brief Method that prefetches storage data to its preferred location. */

    void prefetch() const
    {
        prefetch( mContext );
    }

    /** Will wait for all outstanding asynchronous data transfers. */

    virtual void wait() const = 0;

    /** Allow for additional row compression. */

    void setCompressThreshold( float ratio );

    /** Getter routine for the number of rows. */

    inline IndexType getNumRows() const;

    /** Getter routine for the number of rows. */

    inline IndexType getNumColumns() const;

    virtual void writeAt( std::ostream& stream ) const;

    void resetDiagonalProperty();

    inline bool hasDiagonalProperty() const;

    virtual MatrixStorageFormat getFormat() const = 0;

    /** This method sets storage for the identity matrix
     *
     *  @param[in] n is the size of the square matrix
     */

    virtual void setIdentity( const IndexType n ) = 0;

    /** This method returns the i-th row of the matrix
     *
     * @param[out] row is the destination array that will contain the row
     * @param[in] i is the row index
     *
     */

    virtual void getRow( ContextArray& row, const IndexType i ) const = 0;

    /** This method returns the diagonal of the matrix.
     *
     * @param[out] diagonal is the destination array
     *
     * The values will be stored in diagonal that can be of any type
     * for which implicit conversion is available.
     */

    virtual void getDiagonal( ContextArray& diagonal ) const = 0;

    /** This method sets the diagonal of a matrix storage.
     *
     * Implementation of this routine must be provided by all derived classes.
     */
    virtual void setDiagonal( const ContextArray& diagonal ) = 0;

    virtual void setDiagonal( const Scalar value ) = 0;

    /******************************************************************
     *  Scaling of elements in a matrix                                *
     ******************************************************************/

    /** This method scales all matrix values with a scalar
     *
     * @param[in] value is the source value
     */
    virtual void scale( const Scalar value ) = 0;

    /** This method scales each row of a matrix with a separate value.
     *
     * @param[in] values contains one scale factor for each row.
     *
     * Each row of the matrix is scaled with the corresponding value.
     */
    virtual void scale( const ContextArray& values ) = 0;

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

    virtual void localize( const _MatrixStorage& global, const Distribution& rowDist );

    /** Get the total number of non-zero values in the matrix.
     *
     *  @return total number of non-zero values
     *
     *  An element is considered to be non-zero if its absolute value
     *  is greater equal than mEpsilon. Zero diagonal elements are also
     *  counted if this->hasDiagonalProperty() is given.
     *
     *  This routine does not count zero elements even if they are stored
     *  (e.g. for dense or dia storage data).
     *
     *  This routine gives exactly the same number of elements that will be
     *  the size of the ja and values array when calling buildCSRData.
     *
     *  The default implementation uses the routine buildCSRSizes, but derived
     *  classes should override it if they can do it more efficiently.
     */

    virtual IndexType getNumValues() const;

    inline const LAMAArray<IndexType>& getRowIndexes() const;

    /**
     * @brief Get the number of entries in each row.
     *
     * @param[out] csrIA size array for rows, csrIA.size() == numRows
     *
     * If this->hasDiagonalProperty() is true, diagonal elements are also counted.
     *
     * \code
     *     DenseStorage<double> dense( ... )
     *     LAMAArray<IndexType> ia;
     *     dense.buildCSRSizes( ia );
     * \endcode
     *
     * This routine must be provided by each matrix storage format.
     */

    virtual void buildCSRSizes( LAMAArray<IndexType>& csrIA ) const = 0;

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
     *     LAMAArray<IndexType> ia, ja;
     *     LAMAArray<float> values;
     *     // get non-zero values of dense matrix as float values
     *     dense.buildCSRData( ia, ja, values );
     * \endcode
     *
     * This routine must be provided by each matrix storage format.
     */

    virtual void buildCSRData(
        LAMAArray<IndexType>& csrIA,
        LAMAArray<IndexType>& csrJA,
        ContextArray& csrValues ) const = 0;

    /** Each storage class must provide a routine to set CSR storage data.
     *
     *  @param numRows number of rows
     *  @param numColumns number of columns
     *  @param numValues number of non-zero values
     *  @param csrIA offset array ia for column indexes, size is numRows + 1
     *  @param csrJA are the column indexes of matrix entries, size is numValues
     *  @param csrValues are the values of matrix entries, size is numValues
     */

    virtual void setCSRData(
        const IndexType numRows,
        const IndexType numColumns,
        const IndexType numValues,
        const LAMAArray<IndexType>& csrIA,
        const LAMAArray<IndexType>& csrJA,
        const ContextArray& csrValues ) = 0;

    /** Assign of matrix storage with any format or value type.
     *
     *  Format conversion and type conversion is done implicitly.
     *  Derived clauses might override this method with more efficient solutions.
     *
     *  @param[in] other is the matrix storage to be assigned
     */
    virtual void assign( const _MatrixStorage& other ) = 0;

    /** The opposite routine to assign, for convenience as the other way around is
     *  sometimes more efficient
     */
    virtual void copyTo( _MatrixStorage& other ) const = 0;

    /** Override default assignment operator. */

    _MatrixStorage& operator=( const _MatrixStorage& other );

    /******************************************************************
     *   Help routines (ToDo: -> LAMAArrayUtils ?? )                   *
     ******************************************************************/

    /** Help routines to convert arrays with sizes to offsets and vice versa */

    static void offsets2sizes( LAMAArray<IndexType>& offsets );

    static void offsets2sizes( LAMAArray<IndexType>& sizes, const LAMAArray<IndexType>& offsets );

    /** Utitily to compute the offset array from a sizes array. */

    static IndexType sizes2offsets( LAMAArray<IndexType>& sizes );

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

    /** Each matrix storage must provide a routine that creates a new object
     *  of the same matrix type (same format and same value type).
     *
     *  The implementations in derived classes should use covariant
     *  return types to allow polymorphic creation of objects.
     *
     *  Note: as shared_ptr cannot be used for for covariant return types
     *        this routine returns directly a new created object that should
     *        be wrapped as a shared pointer at calling site.
     */

    virtual _MatrixStorage* clone() const = 0;

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

    /**
     * @brief transformation from matrix type to a csr graph
     *
     * transformation from matrix type to a csr graph,
     * so that (Par)Metis can work with it.
     *
     * @param[out]  adjIA   the ia array of the csr graph
     * @param[out]  adjJA   the ja array of the csr graph
     * @param[out]  vwgt    ToDo
     * @param[in]   comm    Communicator is used to determine number of partitions
     * @param[in]   globalRowIndexes ToDo
     * @param[out]  vtxdist ToDo
     * @since 1.1.0
     */
    virtual void buildCSRGraph(
        IndexType* adjIA,
        IndexType* adjJA,
        IndexType* vwgt,
        CommunicatorPtr comm,
        const IndexType* globalRowIndexes = NULL,
        IndexType* vtxdist = NULL ) const;

protected:

    /** Swaps this with other.
     *
     * swap is protected to avoid accidently wrong swaps of base classes which do not
     * implement their own swap.
     *
     * @param[in,out] other the _MatrixStorage to swap this with
     */
    void swap( _MatrixStorage& other );

    virtual void _assignTranspose( const _MatrixStorage& other );

    virtual void _assign( const _MatrixStorage& other );

    /** Returns the number of bytes needed for the current matrix.
     *
     *  This pure function must be implemented by each derived class. Relevant for the number
     *  of bytes is the current size of the used (LAMA) arrays.
     *
     *  Note: This routine does not tell how many memory is really allocated. Storage data
     *        might be allocated on more than one device. Furthermore, it is possible that
     *        arrays have more memory reserved than needed for its current size.
     */

    virtual size_t getMemoryUsageImpl() const =0;

    IndexType mNumRows; //!< numbers of matrix rows

    IndexType mNumColumns; //!< number of matrix columns

    LAMAArray<IndexType> mRowIndexes; //!< used in case of sparse representation of ia

    float mCompressThreshold; //!< ratio at which compression is done, 0 for never, 1 for always

    bool mDiagonalProperty; //!< if true, diagonal elements are always stored at first position in each row

    SCAI_LOG_DECL_STATIC_LOGGER( logger ) //!< logger for this matrix format

    ContextPtr    mContext;//!< preferred context for the storage

protected:

    /** checkDiagonalProperty checks if the diagonal property of this is full filled.
     *
     * @return true if the diagonal property is fulfilled for the matrix data
     */

    virtual bool checkDiagonalProperty() const = 0;
};

class Distribution;

/** The template class MatrixStorage<ValueType> is the base
 *  class for all matrix storage classes of a given ValueType.
 *
 *  @tparam ValueType is the value type of the matrix values.
 */
template<typename ValueType>
class COMMON_DLL_IMPORTEXPORT MatrixStorage: public _MatrixStorage
{
public:

    typedef ValueType StorageValueType;

    /** Constructor of matrix storage contains dimensions of the matrix. */

    MatrixStorage( const IndexType numRows, const IndexType numColumns );

    /** Overwrite default constructor, same as MatrixStorage(0,0)  */

    MatrixStorage();

    /** Destructor. */

    virtual ~MatrixStorage();

    /** Override _MatrixStorage::create with routine that uses covariant return type. */

    virtual MatrixStorage* clone() const = 0;

    /** Override _MatrixStorage::copy with routine that uses covariant return type. */

    virtual MatrixStorage* copy() const = 0;

    /** Implementation of pure method. */

    virtual common::ScalarType getValueType() const;

    /** Construct a matrix from a dense matrix in row-major order (C-style).
     *  Values of the matrix will be considered as zero if their absolute value is smaller than eps.
     *
     * @param[in] numRows      number of rows
     * @param[in] numColumns   number of columns
     * @param[in] values       the dense matrix values in row-major order (C-style)
     * @param[in] eps          threshold value for non-zero elements
     *
     * Sparse matrix formats will have the diagonal property if numRows == numColums
     */

    template<typename OtherValueType>
    void setRawDenseData(
        const IndexType numRows,
        const IndexType numColumns,
        const OtherValueType values[],
        const ValueType eps = 0.0 );

    void setDenseData(
        const IndexType numRows,
        const IndexType numColumns,
        const ContextArray& values,
        const ValueType eps = 0.0 );

    /**
     * @brief fills matrix storage by csr sparse data.
     *
     * @param[in] numRows    number of rows
     * @param[in] numColumns number of columns
     * @param[in] numValues  the number of stored elements in the matrix
     * @param[in] ia         row pointer of the input csr sparse matrix
     * @param[in] ja         column indexes of the input csr sparse matrix
     * @param[in] values     the data values of the input csr sparse matrix
     */

    template<typename OtherValueType>
    void setRawCSRData(
        const IndexType numRows,
        const IndexType numColumns,
        const IndexType numValues,
        const IndexType* const ia,
        const IndexType* const ja,
        const OtherValueType* const values );

    /** Join local and halo storage back into one storage as needed for NoDistribution.
     *  This matrix storage is used as output matrix.
     *
     *
     * @param[in] localData is matrix storage with local data
     * @param[in] haloData is the matrix storage with halo data
     * @param[in] halo is the communicaton halo, contains also mapping halo to global indexes
     * @param[in] colDist is the distribution that has be used for getting the local data
     * @param[in] keepDiagonalProperty if true routine attempts to keep diagonal property
     *
     * Attention: localData might be aliased with this matrix storage,
     *            while haloData must not
     *
     * Derived classes might use a default implementation that is based on joining
     * CSR data with corresponding conversions.
     */

    virtual void joinHalo(
        const _MatrixStorage& localData,
        const _MatrixStorage& haloData,
        const class Halo& halo,
        const class Distribution& colDist,
        const bool keepDiagonalProperty );

    /** Splitting of matrix storage for halo
     *
     *  @param[out] localData will contain all columns owned by this processor
     *  @param[out] haloData will contain all columns not owned by this processor
     *  @param[out] halo will be the exchange schedule
     *  @param[in]  colDist specifies the distribution of the columns
     *  @param[in]  rowDist optional, localizes the rows before splitting
     */

    virtual void splitHalo(
        MatrixStorage<ValueType>& localData,
        MatrixStorage<ValueType>& haloData,
        class Halo& halo,
        const class Distribution& colDist,
        const class Distribution* rowDist ) const;

    /** Special version of splitHalo where this matrix contains no local
     *  data and where haloData is aliased to this matrix.
     *
     *  This routine is used to translate the non-local column indexes
     *  (required values) to local halo indexes and to set up the
     *  exchange schedule.
     */

    virtual void buildHalo( class Halo& halo, const class Distribution& colDist );

    /** This method build for this matrix the local part of a global matrix.
     *
     *  The row distribution specifies which rows of the global matrix will
     *  be used for this storage locally.
     *
     *  Attention: globalMatrix might be aliased to this storage.
     */

    virtual void localize( const _MatrixStorage& globalData, const class Distribution& rowDist );

    /** This routine builds the full matrix storage for a distributed matrix.
     *
     *  This routine is exactly the inverse routine to the localize routine.
     *  After this operation this matrix storage will contain on each processor
     *  the global matrix data.
     */

    virtual void replicate( const _MatrixStorage& localData, const class Distribution& rowDist );

    /** Get a value of the matrix.
     *
     * @param[in] i is the row index, 0 <= i < mNumRows
     * @param[in] j is the colum index, 0 <= j < mNumRows
     * @throw Exception out-of-range is enabled for ASSERT_LEVEL=DEBUG.
     */

    virtual ValueType getValue( const IndexType i, const IndexType j ) const = 0;

    /**
     *  This method builds CSC sparse data (column sizes, row indexes and data values) for a matrix storage.
     *
     *  For efficiency a derived class might override this method, e.g. for CSR storage.
     */
    virtual void buildCSCData(
        LAMAArray<IndexType>& cscIA,
        LAMAArray<IndexType>& cscJA,
        LAMAArray<ValueType>& cscValues ) const;

    /** Format conversion of matrix storage. A default implementation is provided using CSR data.
     *  Derived clauses might override this method with more efficient solutions.
     */
    virtual void assign( const _MatrixStorage& other );

    /** The opposite routine to assign, for convenience as the other way around is
     *  sometimes more efficient
     */
    virtual void copyTo( _MatrixStorage& other ) const;

    /** Transpose of matrix storage. A default implementation is provided using CSR data.
     *  Derived clauses might override this method with more efficient solutions.
     *
     *  Note: type conversion is not supported here
     */
    virtual void assignTranspose( const MatrixStorage<ValueType>& other );

    /** This methods assigns a redistributed matrix.
     *
     *  @param[in] other is the local part of the matrix data to redistribute
     *  @param[in] redistributor contains source distribution of other and target distribution of this
     *
     */

    virtual void redistribute( const _MatrixStorage& other, const class Redistributor& redistributor );

    /** Special case where other storage is CSR of same type avoids temporary CSR conversion. */

    virtual void redistributeCSR( const CSRStorage<ValueType>& other, const class Redistributor& redistributor );

    /** Build a halo matrix with all rows of required indexes */

    virtual void exchangeHalo(
        const class Halo& halo,
        const MatrixStorage<ValueType>& matrix,
        const class Communicator& comm );

    /** Conversion routine of Compressed Sparse Row data to Compressed Sparse Column.  */

    static void convertCSR2CSC(
        LAMAArray<IndexType>& colIA,
        LAMAArray<IndexType>& colJA,
        LAMAArray<ValueType>& colValues,
        const IndexType numColumns,
        const LAMAArray<IndexType>& rowIA,
        const LAMAArray<IndexType>& rowJA,
        const LAMAArray<ValueType>& rowValues,
        const ContextPtr loc );
    /**
     *  Method that joins rows of another matrix storage
     *
     *   row(i) of out contains all elements of in(k) with rowIndexes[k] == i
     */
    static void joinRows(
        LAMAArray<IndexType>& outIA,
        LAMAArray<IndexType>& outJA,
        LAMAArray<ValueType>& outValues,
        const IndexType numRows,
        const LAMAArray<IndexType>& rowIndexes,
        const LAMAArray<IndexType>& inIA,
        const LAMAArray<IndexType>& inJA,
        const LAMAArray<ValueType>& inValues );

    /**
     * Assignment operator for matrix storage. Using 'assign'
     *
     * @param[in]   other the _MatrixStorage to take a copy from.
     * @return      a reference to this.
     */
    MatrixStorage& operator=( const _MatrixStorage& other );

    /******************************************************************
     *   File I/O for MatrixStorage                                    *
     ******************************************************************/

    virtual void writeToFile(
        const std::string& fileName,
        const File::FileType fileType = File::BINARY,
        const File::DataType dataType = File::INTERNAL,
        const File::IndexDataType indexDataTypeIA = File::INT,
        const File::IndexDataType indexDataTypeJA = File::INT ) const;

    virtual void writeToFile(
        const PartitionId size,
        const PartitionId rank,
        const std::string& fileName,
        const File::FileType fileType = File::BINARY,
        const File::DataType dataType = File::INTERNAL,
        const File::IndexDataType indexDataTypeIA = File::INT,
        const File::IndexDataType indexDataTypeJA = File::INT ) const;

    virtual void readFromFile( const std::string& fileName );

    /******************************************************************
     *   invert                                                        *
     ******************************************************************/

    /** @brief assign this storage with the inverse of another matrix.
     *
     *  @param[in] other is matrix to invert, must be square
     *
     *  Note: other matrix storage can be aliased with this storage.
     */
    virtual void invert( const MatrixStorage<ValueType>& other );

    /******************************************************************
     *   Matrix * ( Vector | Matrix )                                  *
     ******************************************************************/

    /** This method implements result = alpha * thisMatrix * x + beta * y.
     *
     *  Each matrix storage must provide this kind of matrix-vector multiplication.
     *  Default implementation throws exception for non-availability in derived class.
     */

    virtual void matrixTimesVector(
        LAMAArray<ValueType>& result,
        const ValueType alpha,
        const LAMAArray<ValueType>& x,
        const ValueType beta,
        const LAMAArray<ValueType>& y ) const;

    virtual void vectorTimesMatrix(
        LAMAArray<ValueType>& result,
        const ValueType alpha,
        const LAMAArray<ValueType>& x,
        const ValueType beta,
        const LAMAArray<ValueType>& y ) const;

    virtual void matrixTimesVectorN(
        LAMAArray<ValueType>& result,
        const IndexType n,
        const ValueType alpha,
        const LAMAArray<ValueType>& x,
        const ValueType beta,
        const LAMAArray<ValueType>& y ) const;

    /** This method implements result = alpha * thisMatrix * x + beta * y
     *  that is executed asynchronously.
     *
     *  Each matrix storage should provide this kind of matrix-vector multiplication.
     *
     *  A default implementation is provided that does the operation synchronsly
     *  and returns a NoSyncToken.
     */

    virtual tasking::SyncToken* matrixTimesVectorAsync(
        LAMAArray<ValueType>& result,
        const ValueType alpha,
        const LAMAArray<ValueType>& x,
        const ValueType beta,
        const LAMAArray<ValueType>& y ) const;

    virtual tasking::SyncToken* vectorTimesMatrixAsync(
        LAMAArray<ValueType>& result,
        const ValueType alpha,
        const LAMAArray<ValueType>& x,
        const ValueType beta,
        const LAMAArray<ValueType>& y ) const;

    /** Assign this = alpha * a
     *
     *  The default implementation assigns the matrix a and scales it afterwards.
     */

    virtual void matrixTimesScalar( const ValueType alpha, const MatrixStorage<ValueType>& a );

    /** assign alpha * a * b + beta * c */

    virtual void matrixTimesMatrix(
        const ValueType alpha,
        const MatrixStorage<ValueType>& a,
        const MatrixStorage<ValueType>& b,
        const ValueType beta,
        const MatrixStorage<ValueType>& c );

    /** @brief Assign this matrix with alpha * a + beta * b
     *
     *  @param[in] alpha scalar factor for first matrix (storage)
     *  @param[in] a the first matrix (storage)
     *  @param[in] beta scalar factor for second matrix (storage)
     *  @param[in] b the second matrix (storage)
     *
     *  Preconditions:
     *   - a.getNumRows() == b.getNumRows()
     *   - a.getNumColumns() == b.getNumColumns()
     *
     *  \code
     *  result.matrixPlusMatrix( alpha, a, beta, b );
     *  // same as
     *  e.setIdentity( a.getNumColumns() );
     *  result.matrixTimesMatrix( alpha, a, e, beta, b );
     *  \endcode
     *
     *  Location of computation is given by location of first matrix a.
     */

    virtual void matrixPlusMatrix(
        const ValueType alpha,
        const MatrixStorage<ValueType>& a,
        const ValueType beta,
        const MatrixStorage<ValueType>& b );

	/**
     * @brief Returns the L1 norm of this.
     *
     * @return the L1 norm of this.
     *
     * l1Norm computes the sum of the absolute values of this.
     */
    virtual ValueType l1Norm() const = 0;

    /**
     * @brief Returns the L2 norm of this.
     *
     * @return the L2 norm of this.
     *
     * l2Norm computes the sum of the absolute values of this.
     */
    virtual ValueType l2Norm() const = 0;

    /** Get the maximum norm of this matrix
     *
     *  @return maximal absolute value of matrix elements
     */

    virtual ValueType maxNorm() const = 0;

    /** Gets the maximal absolute element-wise difference between two matrices
     *
     *  @param[in] other is the matrix for comparison
     *  @return ( this - other ).maxNorm()
     *
     *  Note: this routine should be more efficient than building the difference matrix
     *        and computing maxNorm of it.
     */

    virtual ValueType maxDiffNorm( const MatrixStorage<ValueType>& other ) const;

    /******************************************************************
     *   Solver methods (e.g. Jacobi )                                 *
     ******************************************************************/

    /** Jacobi iteration step on a local storage
     *
     *  solution = omega * ( rhs + B * oldSolution) * dinv  + ( 1 - omega ) * oldSolution
     *
     *     where B is this storage without diagonal, dinv the inverse diagonal
     *
     *  @param[out] solution is the solution vector to be computed
     *  @param[in]  oldSolution is the old solution vector
     *  @param[in]  rhs is the right hand side of equation system to solve
     *  @param[in]  omega is the scaling factor
     */

    virtual void jacobiIterate(
        LAMAArray<ValueType>& solution,
        const LAMAArray<ValueType>& oldSolution,
        const LAMAArray<ValueType>& rhs,
        const ValueType omega ) const;

    /** Asynchrounous version of jacobiIterate */

    virtual tasking::SyncToken* jacobiIterateAsync(
        LAMAArray<ValueType>& solution,
        const LAMAArray<ValueType>& oldSolution,
        const LAMAArray<ValueType>& rhs,
        const ValueType omega ) const;

    /** Jacobi iteration step on a halo storage.
     *
     *  solution -= omega * ( B(halo) * oldSolution) * dinv
     *
     *  @param[in,out] localSolution is the solution vector that is updated
     *  @param[in]     localStorage is needed to get the diagonal
     *  @param[in]     haloOldSolution is the old solution vector of halo part
     *  @param[in]     omega is the scaling factor.
     *
     *  While local storage must be square and have the diagonal property, this
     *  matrix storage does not have to be square.
     */

    virtual void jacobiIterateHalo(
        LAMAArray<ValueType>& localSolution,
        const MatrixStorage<ValueType>& localStorage,
        const LAMAArray<ValueType>& haloOldSolution,
        const ValueType omega ) const;

    /** Jacobi iteration step on a halo storage.
     *
     *  solution -= omega * ( B(halo) * oldSolution) * dinv
     *
     *  @param[in,out] localSolution is the solution vector that is updated
     *  @param[in]     localDiagonal pointer to the diagonal of local storage
     *  @param[in]     haloOldSolution is the old solution vector of halo part
     *  @param[in]     omega is the scaling factor.
     *
     *  @since 1.1.0
     */

    virtual void jacobiIterateHalo(
        LAMAArray<ValueType>& localSolution,
        const LAMAArray<ValueType>& localDiagonal,
        const LAMAArray<ValueType>& haloOldSolution,
        const ValueType omega ) const;

    // Note: Asynchronous version of jacobiIterateHalo not supported

    using _MatrixStorage::getContext;
    using _MatrixStorage::getContextPtr;

    // Use this method to change epsiolon temporarily

    void swapEpsilon( ValueType& epsilon ) const
    {
        std::swap( epsilon, mEpsilon );
    }

    /**
     * @brief Returns whether the matrix storage is symmetric or not.
     *
     * @return a boolean pointing out whether the matrix storage is symmetric or not.
     */
    bool checkSymmetry() const;

protected:

    /** The value mEpsilon is an individual value for each matrix storage that
     *  specifies a threshold when a matrix values can be considered as zero.
     *  It is used internally especially when setting dense data or when
     *  the storage data is compressed.
     */

    mutable ValueType mEpsilon;

    /** Swap member variables of base class MatrixStorage<ValueType>
     *
     * @param[in,out] other the MatrixStorage to swap this with
     *
     * swap is protected to avoid accidently wrong swaps of derived classes which do not
     * implement their own swap or where storages of different format or types are involved.
     *
     * \code
     * CSRStorage<float> csr;
     * CSRStorage<double> csr1;
     * ELLStorage<float> ell;
     * csr.swap( ell ); // NOT ALLOWED, different format
     * csr.swap( csr1 ); // NOT ALLOWED, different type
     * \endcode
     */
    void swap( MatrixStorage<ValueType>& other );
};

/* ------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------- */

ContextPtr _MatrixStorage::getContextPtr() const
{
    return mContext;
}

const Context& _MatrixStorage::getContext() const
{
    return *mContext;
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

const LAMAArray<IndexType>& _MatrixStorage::getRowIndexes() const
{
    return mRowIndexes;
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
template<typename OtherValueType>
void MatrixStorage<ValueType>::setRawCSRData(
    const IndexType numRows,
    const IndexType numColumns,
    const IndexType numValues,
    const IndexType* const ia,
    const IndexType* const ja,
    const OtherValueType* const values )
{
    // wrap the pointer data into LAMA arrays ( without copies )

    LAMAArrayRef<IndexType> csrIA( numRows + 1, ia );
    LAMAArrayRef<IndexType> csrJA( numValues, ja );
    LAMAArrayRef<OtherValueType> csrValues( numValues, values );

    // now set the data on the context of this storage via virtual method

    setCSRData( numRows, numColumns, numValues, csrIA, csrJA, csrValues );
}

} /* end namespace lama */

} /* end namespace scai */

/* ------------------------------------------------------------------------- */

COMMON_DLL_IMPORTEXPORT std::ostream& operator<<( std::ostream& stream, const lama::MatrixStorageFormat storageFormat );

