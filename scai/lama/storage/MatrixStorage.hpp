/**
 * @file MatrixStorage.hpp
 *
 * @license
 * Copyright (c) 2009-2018
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the SCAI framework LAMA.
 *
 * LAMA is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * LAMA is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
 * more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with LAMA. If not, see <http://www.gnu.org/licenses/>.
 * @endlicense
 *
 * @brief Definition of a common base class for all matrix storage formats.
 * @author Thomas Brandes
 * @date 27.04.2011
 */

#pragma once

// local library
#include <scai/lama/storage/_MatrixStorage.hpp>

#include <scai/lama/freeFunction.hpp>

#include <scai/lama/io/FileIO.hpp>
#include <scai/lama/storage/Format.hpp>

#include <scai/dmemo/Communicator.hpp>
#include <scai/hmemo.hpp>

// internal scai libraries
#include <scai/common/BinaryOp.hpp>
#include <scai/common/UnaryOp.hpp>
#include <scai/common/MatrixOp.hpp>
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

template<typename ValueType> class CSRStorage;

/** The template class MatrixStorage<ValueType> is the base
 *  class for all matrix storage classes of a given ValueType.
 *
 *  It is derived from the non-templace class _MatrixStorage and
 *  does not add any new member variable.
 *
 *  @tparam ValueType is the value type of the matrix values.
 */
template<typename ValueType>
class COMMON_DLL_IMPORTEXPORT MatrixStorage: public _MatrixStorage
{
public:

    typedef ValueType StorageValueType;

    /** Destructor. */

    virtual ~MatrixStorage();

    /** Override _MatrixStorage::copy with routine that uses covariant return type. */

    virtual MatrixStorage* copy() const = 0;

    /** Override _MatrixStorage::newMatrixStorage with routine that uses covariant return type. */

    virtual MatrixStorage* newMatrixStorage( const IndexType numRows, const IndexType numColumns ) const = 0;

    MatrixStorage* newMatrixStorage() const
    {
        return newMatrixStorage( getNumRows(), getNumColumns() );
    }

    /** Implementation of pure method. */

    virtual common::ScalarType getValueType() const;

    /** Construct a matrix stroage from raw data in row-major order (C-style).
     *  Values of the matrix will be considered as zero if their value is exact zero.
     *
     * @param[in] numRows      number of rows
     * @param[in] numColumns   number of columns
     * @param[in] values       the dense matrix values in row-major order (C-style)
     */
    template<typename OtherValueType>
    void setRawDenseData(
        const IndexType numRows,
        const IndexType numColumns,
        const OtherValueType values[] );

    /** This method scales all matrix values with a scalar
     *
     * @param[in] value is the source value
     */
    virtual void scale( const ValueType value ) = 0;

    virtual void conj() = 0;

    /******************************************************************
     *  Scaling of elements in a matrix                                *
     ******************************************************************/

    /** This method scales each row of a matrix with a separate value.
     *
     * @param[in] values contains one scale factor for each row.
     *
     * Each row of the matrix is scaled with the corresponding value.
     */
    virtual void scaleRows( const hmemo::HArray<ValueType>& values ) = 0;

    virtual void reduce( 
        hmemo::HArray<ValueType>& array, 
        const IndexType dim, 
        const common::BinaryOp reduceOp, 
        const common::UnaryOp elemOp );

    /**
     * @brief fills any matrix storage by raw csr sparse data.
     *
     * @param[in] numRows    number of rows,
     * @param[in] numColumns number of columns
     * @param[in] numValues  the number of stored elements in the matrix
     * @param[in] ia         csr offset array, size is numRows + 1
     * @param[in] ja         column indexes of the input csr sparse matrix, size is numValues
     * @param[in] values     the data values of the input csr sparse matrix, size is numValues
     */
    template<typename OtherValueType>
    void setRawCSRData(
        const IndexType numRows,
        const IndexType numColumns,
        const IndexType numValues,
        const IndexType* const ia,
        const IndexType* const ja,
        const OtherValueType* const values );

    /* =========================================================== */
    /*    set/get for rows + columns                               */
    /* =========================================================== */

    /** This method returns the i-th row of the matrix
     *
     * @param[out] row is the destination array that will contain the row
     * @param[in] i is the row index
     */
    virtual void getRow( hmemo::HArray<ValueType>& row, const IndexType i ) const = 0;

    /** This method updates the i-th row of the matrix
     *
     * @param[in] row is the source array that will contain the update values
     * @param[in] i is the row index
     * @param[in] op specifies how the elements should e combined
     *
     * Note: in sparse matrices only existing entries are updated.
     */
    virtual void setRow( const hmemo::HArray<ValueType>& row, const IndexType i, const common::BinaryOp op ) = 0;

    /** This method returns the i-th row of the matrix as a sparse array, i.e. indexes + values
     *
     * @param[out] jA are column positions with non-zero values
     * @param[out] values are the values at the corresponding column positions.
     * @param[in] i is the row index
     *
     */
    virtual void getSparseRow( hmemo::HArray<IndexType>& jA, hmemo::HArray<ValueType>& values, const IndexType i ) const = 0;

    /** This method returns the j-th col of the matrix as a sparse array, i.e. indexes + values
     *
     * @param[out] iA are row positions with non-zero values
     * @param[out] values are the values at the corresponding row positions.
     * @param[in] j is the index for queried column
     *
     */
    virtual void getSparseColumn( hmemo::HArray<IndexType>& iA, hmemo::HArray<ValueType>& values, const IndexType j ) const = 0;

    virtual void setColumn( const hmemo::HArray<ValueType>& column, const IndexType j,
                            common::BinaryOp op            ) = 0;

    /** This method returns the j-th column of the matrix
     *
     * @param[out] column is the destination array that will contain the column
     * @param[in] j is the row index
     *
     */
    virtual void getColumn( hmemo::HArray<ValueType>& column, const IndexType j ) const = 0;

    /* =========================================================== */
    /*    set/get for diagonals                                    */
    /* =========================================================== */

    /**
     * @brief Assign this matrix a diagonal matrix specified by an array with the diagonal elements.
     *
     * @param[in] diagonal contains the values for the diagonal
     *
     * \code
     *     MatrixStorage<ValueType>& m = ...
     *     // m.setIdentity( n ) can also be written as follows
     *     m.assignDiagonal( HArray<ValueType>( n, ValueType( 1 ) ) );  
     * \endcode
     */
    virtual void assignDiagonal( const hmemo::HArray<ValueType>& diagonal ) = 0;

    /** This method sets the diagonal of a matrix storage.
     *
     * Implementation of this routine must be provided by all derived classes.
     * It might throw an exception if the diagonal property is not given 
     * (i.e. the sparse pattern has not all entries for the diagonal elements).
     */
    virtual void setDiagonalV( const hmemo::HArray<ValueType>& diagonal ) = 0;

    /**
     * This method sets the diagonal of a matrix storage with the same single value.
     *
     * Implementation of this routine must be provided by all derived classes.
     * It might throw an exception if the diagonal property is not given 
     * (i.e. the sparse pattern has not all entries for the diagonal elements).
     */
    virtual void setDiagonal( const ValueType value ) = 0;

    /** This method returns the diagonal of this matrix storage.
     *
     * @param[out] diagonal will contain the diagonal of this matrix.
     *
     * The size of the output array will be min( nRows, nCols ).
     *
     * This pure method must be implemented by all storage classes. 
     */
    virtual void getDiagonal( hmemo::HArray<ValueType>& diagonal ) const = 0;

    /** Join local and halo storage back into one storage as needed for NoDistribution.
     *  This matrix storage is used as output matrix.
     *
     * @param[in] localData is matrix storage with local data
     * @param[in] haloData is the matrix storage with halo data
     * @param[in] halo is the communicaton halo, contains also mapping halo to global indexes
     * @param[in] colDist is the distribution that has be used for getting the local data
     *
     * Attention: localData might be aliased with this matrix storage,
     *            while haloData must not
     *
     * Derived classes might use a default implementation that is based on joining
     * CSR data with corresponding conversions.
     */
    virtual void joinHalo(
        const MatrixStorage<ValueType>& localData,
        const MatrixStorage<ValueType>& haloData,
        const dmemo:: Halo& halo,
        const dmemo::Distribution& colDist );

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
        dmemo::Halo& halo,
        const dmemo::Distribution& colDist,
        const dmemo::Distribution* rowDist ) const;

    /** Special version of splitHalo where this matrix contains no local
     *  data and where haloData is aliased to this matrix.
     *
     *  This routine is used to translate the non-local column indexes
     *  (required values) to local halo indexes and to set up the
     *  exchange schedule.
     */

    virtual void buildHalo( dmemo::Halo& halo, const dmemo::Distribution& colDist );

    /** This method translates the halo column indexes back to global indexes 
     * 
     *  @param halo provides the method for translation of halo to global indexes
     *  @param globalNumColumns will be the new number of columns for this matrix
     * 
     *  The default implementation uses CSR storage to globalize the indexes.
     */
    virtual void globalizeHaloIndexes( const dmemo::Halo& halo, const IndexType globalNumColumns );

    /**
     * @brief This method removes all zero elements of a sparse storage, i.e. only entries whose absolute
     *        value is greater than eps are considered to be non-zero.
     *
     * @param[in] eps  is the threshold when a values is to be considered as zero
     *
     * The default implementation uses temporary CSR storage to compress it.
     *
     * \code
     *    auto diffStorage = eval<CSRStorage<ValueType>>( storage1 - storage2 );
     *    diffStorage.compress( 0.0001 );
     * \endcode
     */
    virtual void compress( const RealType<ValueType> eps = 0 );

    /** This method build for this matrix the local part of a global matrix.
     *
     *  The row distribution specifies which rows of the global matrix will
     *  be used for this storage locally.
     *
     *  Attention: globalMatrix might be aliased to this storage.
     */

    virtual void localize( const _MatrixStorage& globalData, const dmemo::Distribution& rowDist );

    /** This routine builds the full matrix storage for a distributed matrix.
     *
     *  This routine is exactly the inverse routine to the localize routine.
     *  After this operation this matrix storage will contain on each processor
     *  the global matrix data.
     */

    virtual void replicate( const _MatrixStorage& localData, const dmemo::Distribution& rowDist );

    /** Get a value of the matrix.
     *
     * @param[in] i is the row index, 0 <= i < getNumRows()
     * @param[in] j is the colum index, 0 <= j < getNumColumns()
     * @throw Exception out-of-range is enabled for ASSERT_LEVEL=DEBUG.
     */

    virtual ValueType getValue( const IndexType i, const IndexType j ) const = 0;

    /** Set/update an existing value of the matrix.
     *
     *  @param[in] i is the row index, 0 <= i < getNumRows()
     *  @param[in] j is the col index, 0 <= j < getNumColumns()
     *  @param[in] val is the value to be set
     *  @param[in] op specifies how to combine old and new value
     *
     *  @throw Exception if value is non-zero and sparse pattern does not contain element
     */
    virtual void setValue( const IndexType i,
                           const IndexType j,
                           const ValueType val,
                           const common::BinaryOp op = common::BinaryOp::COPY ) = 0;

    /**
     *  This method builds CSC sparse data (column sizes, row indexes and data values) for a matrix storage.
     *
     *  For efficiency a derived class might override this method, e.g. for CSR storage.
     */
    virtual void buildCSCData(
        hmemo::HArray<IndexType>& cscIA,
        hmemo::HArray<IndexType>& cscJA,
        hmemo::HArray<ValueType>& cscValues ) const;

    /** Format conversion of matrix storage. A default implementation is provided using CSR data.
     *  Derived clauses might override this method with more efficient solutions.
     */
    virtual void assignDummy( const _MatrixStorage& other );

    /** The opposite routine to assign, for convenience as the other way around is
     *  sometimes more efficient
     */
    virtual void copyTo( _MatrixStorage& other ) const;

    /**
     *  Implementation of _MatrixStorage::copyBlockTo
     */
    virtual void copyBlockTo( _MatrixStorage& other, const IndexType first, const IndexType n ) const;

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

    virtual void redistribute( const _MatrixStorage& other, const dmemo::Redistributor& redistributor );

    /** Special case where other storage is CSR of same type avoids temporary CSR conversion. */

    virtual void redistributeCSR( const CSRStorage<ValueType>& other, const dmemo::Redistributor& redistributor );

    /** Build a halo matrix with all rows of required indexes */

    virtual void exchangeHalo(
        const dmemo::Halo& halo,
        const MatrixStorage<ValueType>& matrix,
        const dmemo::Communicator& comm );

    /**
     *  Method that joins rows of another matrix storage
     *
     *   row(i) of out contains all elements of in(k) with rowIndexes[k] == i
     */
    static void joinRows(
        hmemo::HArray<IndexType>& outIA,
        hmemo::HArray<IndexType>& outJA,
        hmemo::HArray<ValueType>& outValues,
        const IndexType numRows,
        const hmemo::HArray<IndexType>& rowIndexes,
        const hmemo::HArray<IndexType>& inIA,
        const hmemo::HArray<IndexType>& inJA,
        const hmemo::HArray<ValueType>& inValues );

    /**
     * Assignment operator for matrix storage. Using 'assign'
     *
     * @param[in]   other the _MatrixStorage to take a copy from.
     * @return      a reference to this.
     */
    MatrixStorage& operator=( const _MatrixStorage& other );

    MatrixStorage& operator=( const MatrixStorage<ValueType>& other );

    /** Implementation of pure routine _MatrixStorage::writeToFile */

    virtual void writeToFile(
        const std::string& fileName,
        const std::string& fileType = "",
        const common::ScalarType dataType = common::ScalarType::UNKNOWN,
        const common::ScalarType indexType = common::ScalarType::UNKNOWN,
        const FileIO::FileMode fileMode = FileIO::DEFAULT_MODE  ) const;

    virtual void writeToFile(
        const PartitionId size,
        const PartitionId rank,
        const std::string& fileName,
        const std::string& fileType,
        const common::ScalarType dataType = common::ScalarType::UNKNOWN,
        const common::ScalarType indexType = common::ScalarType::UNKNOWN,
        const FileIO::FileMode fileMode = FileIO::DEFAULT_MODE  ) const;

    virtual void readFromFile( const std::string& fileName, const IndexType firstRow = 0, const IndexType nRows = invalidIndex );

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
     *   _Matrix * ( Vector | _Matrix )                                  *
     ******************************************************************/

    /** This method implements result = alpha * thisMatrix * x + beta * y.
     *
     *  Each matrix storage must provide this kind of matrix-vector multiplication.
     */

    virtual void matrixTimesVector(
        hmemo::HArray<ValueType>& result,
        const ValueType alpha,
        const hmemo::HArray<ValueType>& x,
        const ValueType beta,
        const hmemo::HArray<ValueType>& y,
        const common::MatrixOp op ) const = 0;

    virtual void matrixTimesVectorN(
        hmemo::HArray<ValueType>& result,
        const IndexType n,
        const ValueType alpha,
        const hmemo::HArray<ValueType>& x,
        const ValueType beta,
        const hmemo::HArray<ValueType>& y ) const;

    /** This method implements result = alpha * thisMatrix * x + beta * y
     *  that is executed asynchronously.
     *
     *  Each matrix storage should provide this kind of matrix-vector multiplication.
     *
     *  A default implementation is provided that does the operation synchronsly
     *  and returns a NoSyncToken.
     */

    virtual tasking::SyncToken* matrixTimesVectorAsync(
        hmemo::HArray<ValueType>& result,
        const ValueType alpha,
        const hmemo::HArray<ValueType>& x,
        const ValueType beta,
        const hmemo::HArray<ValueType>& y,
        const common::MatrixOp op ) const;

    /** Assign this = alpha * a
     *
     *  The default implementation assigns the matrix a and scales it afterwards.
     */

    virtual void matrixTimesScalar( const ValueType alpha, const MatrixStorage<ValueType>& a );

    /** this = alpha * a * b + beta * c 
     *
     *  This class provides a default implementation where it calls the routine for CSR storage
     *  and converts the result afterwards.
     */
    virtual void matrixTimesMatrix(
        const ValueType alpha,
        const MatrixStorage<ValueType>& a,
        const MatrixStorage<ValueType>& b,
        const ValueType beta,
        const MatrixStorage<ValueType>& c );
 
    /** Apply binary operation elementwise for matrix storage elements. */

    virtual void binaryOp(
        const MatrixStorage<ValueType>& a,
        const common::BinaryOp op,
        const MatrixStorage<ValueType>& b );
       
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
    virtual RealType<ValueType> l1Norm() const = 0;

    /**
     * @brief Returns the L2 norm of this.
     *
     * @return the L2 norm of this.
     *
     * l2Norm computes the sum of the absolute values of this.
     */
    virtual RealType<ValueType> l2Norm() const = 0;

    /** Get the maximum norm of this matrix
     *
     *  @return maximal absolute value of matrix elements
     */

    virtual RealType<ValueType> maxNorm() const = 0;

    /** Gets the maximal absolute element-wise difference between two matrices
     *
     *  @param[in] other is the matrix for comparison
     *  @return ( this - other ).maxNorm()
     *
     *  Note: this routine should be more efficient than building the difference matrix
     *        and computing maxNorm of it.
     */

    virtual RealType<ValueType> maxDiffNorm( const MatrixStorage<ValueType>& other ) const;

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
        hmemo::HArray<ValueType>& solution,
        const hmemo::HArray<ValueType>& oldSolution,
        const hmemo::HArray<ValueType>& rhs,
        const ValueType omega ) const;

    /** Asynchrounous version of jacobiIterate */

    virtual tasking::SyncToken* jacobiIterateAsync(
        hmemo::HArray<ValueType>& solution,
        const hmemo::HArray<ValueType>& oldSolution,
        const hmemo::HArray<ValueType>& rhs,
        const ValueType omega ) const;

    /** Jacobi iteration step on a halo storage.
     *
     *  solution -= omega * ( B(halo) * oldSolution) ./ localDiagonal
     *
     *  @param[in,out] localSolution is the solution vector that is updated
     *  @param[in]     localDiagonal pointer to the diagonal of local storage
     *  @param[in]     haloOldSolution is the old solution vector of halo part
     *  @param[in]     omega is the scaling factor.
     */
    virtual void jacobiIterateHalo(
        hmemo::HArray<ValueType>& localSolution,
        const hmemo::HArray<ValueType>& localDiagonal,
        const hmemo::HArray<ValueType>& haloOldSolution,
        const ValueType omega ) const;

    // Note: Asynchronous version of jacobiIterateHalo not supported

    using _MatrixStorage::getContextPtr;

    /**
     * @brief Returns whether the matrix storage is symmetric or not.
     *
     * @return a boolean pointing out whether the matrix storage is symmetric or not.
     */
    virtual bool checkSymmetry() const;

    /* ========================================================================= */
    /*       Filling matrix storage with assembled COO data                      */
    /* ========================================================================= */

    /**
     *  @brief Add assembled data in coordinate format (COO) to the matrix storage
     *
     *  @param[in] ia     array with row positions
     *  @param[in] ja     array with col positions, same size as ia
     *  @param[in] values same size as ia, ja, values[k] is value for ( ia[k], ja[k] )
     *  @param[in] op     specifies how to combine with existing entries.
     *
     *  A default implementation is provided here, derived classes might override it.
     */
    virtual void fillCOO( 
        hmemo::HArray<IndexType> ia, 
        hmemo::HArray<IndexType> ja, 
        hmemo::HArray<ValueType> values,
        const common::BinaryOp op = common::BinaryOp::COPY );
               
protected:

    /** Constructor of matrix storage contains dimensions of the matrix. */

    MatrixStorage( const IndexType numRows, const IndexType numColumns, hmemo::ContextPtr ctx );

    /** Move constructor that can be used by derived classes */

    MatrixStorage( MatrixStorage<ValueType>&& other ) noexcept;

    /** Copy constructor that can be used by derived classes. */

    MatrixStorage( const MatrixStorage<ValueType>& other );

    /** Move assignment operator, explicity deleted as there is no virtual move semantic.
     */
    MatrixStorage& operator=( MatrixStorage<ValueType>&& other ) = delete;

    /** Method that provides move semantic for this base class.
     
     *  Note: this method only moves member variables of this base class but
     *        not those of the derived classes.
     */
    void moveImpl( MatrixStorage<ValueType>&& other );

    /** Common routine that checks sizes of matrixTimesVector call */

    void gemvCheck(
        const ValueType alpha,
        const hmemo::HArray<ValueType>& x,
        const ValueType beta,
        const hmemo::HArray<ValueType>& y,
        const common::MatrixOp op ) const;

public:

    static MatrixStorage<ValueType>* create( const MatrixStorageCreateKeyType key );

};

/* ------------------------------------------------------------------------- */
/* Implementation of inline methods                                          */
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
    hmemo::HArrayRef<IndexType> csrIA( numRows + 1, ia );
    hmemo::HArrayRef<IndexType> csrJA( numValues, ja );
    hmemo::HArrayRef<OtherValueType> csrValues( numValues, values );
    // now set the data on the context of this storage via virtual method
    setCSRData( numRows, numColumns, csrIA, csrJA, csrValues );
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
MatrixStorage<ValueType>* MatrixStorage<ValueType>::create( const MatrixStorageCreateKeyType key )
{
    return static_cast<MatrixStorage<ValueType>* >( _MatrixStorage::create( key ) );
}

} /* end namespace lama */

} /* end namespace scai */
