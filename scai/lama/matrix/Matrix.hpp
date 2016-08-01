/**
 * @file Matrix.hpp
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
 * @brief Abstract base class for all matrices supported by LAMA.
 * @author Jiri Kraus, Thomas Brandes
 * @date 22.02.2011
 */
#pragma once

// for dll_import
#include <scai/common/config.hpp>

// base classes
#include <scai/dmemo/Distributed.hpp>

// local library
#include <scai/lama/Scalar.hpp>
#include <scai/lama/Vector.hpp>
#include <scai/lama/io/FileIO.hpp>

#include <scai/dmemo/Distribution.hpp>
#include <scai/dmemo/NoDistribution.hpp>

#include <scai/lama/expression/Expression.hpp>
#include <scai/lama/storage/MatrixStorage.hpp>

// internal scai libraries
#include <scai/hmemo.hpp>

#include <scai/logging.hpp>

#include <scai/common/Factory.hpp>
#include <scai/common/ScalarType.hpp>
#include <scai/common/SCAITypes.hpp>

namespace scai
{

namespace lama
{

/** Pointer class for a matrix, always use of a shared pointer. */

typedef common::shared_ptr<class Matrix> MatrixPtr;


/** Key type used for the Matrix factory.
 *
 *  Note: own struct instead of std::pair to allow definition of operator <<
 */

typedef MatrixStorageCreateKeyType MatrixCreateKeyType;

/**
 * @brief The class Matrix is a abstract type that represents a distributed 2D real or complex matrix.
 *
 * Matrix is one of the LAMA Base Types and should be used in all situations where it is not necessary to access a
 * single element or to create a new Matrix.
 */
class COMMON_DLL_IMPORTEXPORT Matrix:

    public common::Factory<MatrixCreateKeyType, Matrix*>,
    public dmemo::Distributed,
    public Format

{

public:

    /** @brief More convenient use of the create routine of factory that avoids use of CreateKeyType.
     */
    static Matrix* getMatrix( const MatrixStorageFormat format, const common::scalar::ScalarType valueType );

    /**
     * @brief ExpressionMemberType is the type that is used the template Expression to store a Vector.
     */
    typedef const Matrix& ExpressionMemberType;

    /**
     * @brief Destructor, releases all allocated resources.
     */
    virtual ~Matrix();

    /** Override Distributed::buildCSRGraph */

    virtual void buildCSRGraph( IndexType ia[], IndexType ja[], IndexType vwgt[], const IndexType* globalRowIndexes ) const;

    /** Override Distributed::getCSRGraphSize */

    virtual IndexType getCSRGraphSize() const;

    /**
     * @brief write the matrix to an output file
     *
     * @param[in] fileName is the name of the output file (suffix must be added according to the file type)
     * @param[in] fileType format of the output file (SAMG, MatrixMarket), default is to decide by suffix
     * @param[in] valuesType representation type for output values, default is same type as matrix values
     * @param[in] indexType representation type for col/row index values
     * @param[in] fileMode can be used to forche BINARY or FORMATTED output
     */

    void writeToFile(
        const std::string& fileName,
        const std::string& fileType = "",
        const common::scalar::ScalarType dataType = common::scalar::UNKNOWN,
        const common::scalar::ScalarType indexType = common::scalar::UNKNOWN,
        const FileIO::FileMode fileMode = FileIO::DEFAULT_MODE  ) const;

    /**
     * @brief Checks for a given matrix whether the content of its data is sound.
     *
     * @return false if any of the internal data structures is not okay
     *
     * This method returns the same value on all processors.
     *
     * If any inconsistency has been found an error message should be logged, but it should
     * not throw an exception. This might be done by the caller of this routine to avoid
     * working with inconsistent matrices.
     *
     * \code
     * SCAI_ASSERT_DEBUG( a.isConsistent(), a << ": is invalid matrix after reading" )
     * \endcode
     */

    virtual bool isConsistent() const = 0;

    /**
     *  @brief Virtual method that delivers the class name to which a matrix belongs.
     */
    virtual const char* getTypeName() const = 0;

    /**
     * @brief Clears the full matrix, resets global and local sizes to 0.
     *
     * \code
     *     CSRSparseMatrix<double> a ( ... )
     *     a = CSRSparseMatrix<double> ();     \\ will free all arrays
     *     a.clear();                          \\ same functionality, clears involved arrays
     *
     * \endcode
     *
     * This routine will not free any allocated memory. It is especially helpful to invalidate
     * all data before it is defined again.
     */
    virtual void clear() = 0;

    /**
     * @brief Same as clear but here all data is freed.
     */
    virtual void purge() = 0;

    /** @brief Reallocates this matrix to a replicated zero-matrix of the given shape.
     *
     * @param[in] numRows      number of rows, must be non-negative.
     * @param[in] numColumns   number of columns, must be non-negative.
     *
     *  @remark The call allocate( 0, 0 ) implies a clear on all arrays used internally
     *          for the presentation of the matrix data.
     */
    virtual void allocate( const IndexType numRows, const IndexType numColumns ) = 0;

    /** @brief Reallocates this matrix to a distributed zero-matrix by the given distributions.
     *
     *  @param[in] rowDistribution is row distribution, number of rows given by getGlobalSize()
     *  @param[in] colDistribution is col distribution, number of columns given by getGlobalSize()
     *
     *  The distributions specify the new global and the new local sizes for a distributed matrix.
     */
    virtual void allocate( dmemo::DistributionPtr rowDistribution, dmemo::DistributionPtr colDistribution ) = 0;

    /**
     * @brief Operator that sets the matrix to the identity matrix.
     *
     * \code
     * void sub( ..., Matrix& a, ... )
     * ...
     * SCAI_ASSERT_EQUAL_DEBUG( a.getNumRows(), a.getNumColumns() )
     * a.setIdentity( a.getRowDistribution() );
     * \endcode
     */

    /** Set matrix to a identity square matrix with same row and column distribution. */

    virtual void setIdentity( dmemo::DistributionPtr distribution ) = 0;

    /** Set matrix to a (replicated) identity matrix with same row and column distribution. 
     *
     *  \code
     *    m.setIdentitiy( n ) ->  m.setIdentity( DistributionPtr( new NoDistribution( n ) ) );
     *  \code
     */

    void setIdentity( const IndexType n );

    /**
     * This method sets a matrix by reading its values from one or multiple files.
     *
     * @param[in] filename      the filename to read from
     * @param[in] rowDist       optional, if set it is the distribution of the matrix 
     *
     *   \code
     *      CSRSparseMatrix<double> matrix;
     *      matrix.readFromFile( "matrix.mtx" )                    ! matrix only on processor 0
     *      matrix.readFromFile( "matrix_%r.mtx" )                 ! general block distributed matrix, each processor reads it own file
     *      matrix.readFromFile( "matrix.mtx", rowDist )           ! each processor gets its local part of the matrix in one file
     *      matrix.readFromFile( "matrix_%r.mtx", rowDist )        ! read a partitioned matrix with the given distribution
     *   \endcode
     */
    void readFromFile( const std::string& fileName, dmemo::DistributionPtr rowDist = dmemo::DistributionPtr() );

    /**
     *  This method sets a matrix a reading its values from one or multiple files and also the distribution from a file
     *
     * @param[in] matrixFileName the single or partitioned filename to read from
     * @param[in] distributionFileName the single or partitioned filename with the row distribution of the matrix
     *
     *   \code
     *      CSRSparseMatrix<double> matrix;
     *      matrix.readFromFile( "matrix.mtx", "owners.mtx" )
     *      matrix.readFromFile( "matrix_%r.mtx", "owners.mtx" )
     *      matrix.readFromFile( "matrix.mtx", "rows%r.mtx" )
     *      matrix.readFromFile( "matrix_%r.mtx", "rows%r.mtx" )
     *   \endcode
     */
    void readFromFile( const std::string& matrixFileName, const std::string& distributionFileName );

    /** This method resorts column indexes in such a way that the diagonal element is always the 
     *  first one in a row. 
     *
     *  This method throws an exception if row and column distribution are not equal. Furhtermore
     *  it throws an exception, if a diagonal element is zero, i.e. there is no entry for the diagonal
     *  element in a sparse format.
     */
    void setDiagonalProperty();

    /** This method sets a matrix with the values owned by this partition in dense format
     *
     *  @param[in] rowDist distributon of rows for the matrix
     *  @param[in] colDist distributon of columns for the matrix
     *  @param[in] values contains all values of the owned rows in row-major order (C-style)
     *  @param[in] eps    threshold value for non-zero elements
     *
     *  Note: only the row distribution decides which data is owned by this processor
     *
     *  The following must be valid: values.size() == rowDist->getLocalSize() * colDist->getGlobalSize()
     */
    virtual void setDenseData( dmemo::DistributionPtr rowDist, dmemo::DistributionPtr colDist, const hmemo::_HArray& values, Scalar eps =
                                   Scalar( 0 ) ) = 0;

    /** This method set a matrix with the values owned by this partition in CSR format
     *
     *  @param[in] rowDist distributon of rows for the matrix
     *  @param[in] colDist distributon of columns for the matrix
     *  @param[in] numValues number of non-zero values
     *  @param[in] ia     is the offset array for number of elements in local rows
     *  @param[in] ja     contains the (global) column indexes
     *  @param[in] values contains the matrix values for the entries specified by ja
     *
     *  Note: only the row distribution decides which data is owned by this processor
     *
     *  - ja.size() == values.size() must be valid, stands for the number of non-zero values of this partition
     *  - ia.size() == rowDistribution.getLocalSize() + 1 must be valid
     */

    virtual void setCSRData(
        dmemo::DistributionPtr rowDist,
        dmemo::DistributionPtr colDist,
        const IndexType numValues,
        const hmemo::HArray<IndexType>& ia,
        const hmemo::HArray<IndexType>& ja,
        const hmemo::_HArray& values ) = 0;

    /** This method set a matrix with the values owned by this partition in DIA format
     *
     *  @param[in] rowDist distributon of rows for the matrix
     *  @param[in] colDist distributon of columns for the matrix
     *  @param[in] numDiagonals number of stored diagonals
     *  @param[in] offsets offsets of the stored diagonals to the main diagonal
     *  @param[in] values contains the matrix values for each diagonal
     *
     *  Note: only the row distribution decides which data is owned by this processor
     *
     *  - numDiagonals == offset.size() must be valid, stands for the number stored diagonals
     */

    virtual void setDIAData(
        dmemo::DistributionPtr rowDist,
        dmemo::DistributionPtr colDist,
        const IndexType numDiagonals,
        const hmemo::HArray<IndexType>& offsets,
        const hmemo::_HArray& values ) = 0;

    /** This method sets raw dense data in the same way as setDenseData but with raw value array */

    template<typename ValueType>
    void setRawDenseData(
        dmemo::DistributionPtr rowDist,
        dmemo::DistributionPtr colDist,
        const ValueType* values,
        const ValueType eps = 0 )
    {
        const IndexType n = rowDist->getLocalSize();
        const IndexType m = colDist->getGlobalSize();
        // use of HArrayRef instead of HArray avoids additional copying of values
        const hmemo::HArrayRef<ValueType> valueArray( n * m, values );
        setDenseData( rowDist, colDist, valueArray, Scalar( eps ) );
    }

    /** This method sets raw CSR data in the same way as setCSRData but with raw value array */

    template<typename ValueType>
    void setRawCSRData(
        dmemo::DistributionPtr rowDist,
        dmemo::DistributionPtr colDist,
        const IndexType numValues,
        const IndexType* ia,
        const IndexType* ja,
        const ValueType* values )
    {
        const IndexType n = rowDist->getLocalSize();
        // use of HArrayRef instead of HArray avoids additional copying of values
        const hmemo::HArrayRef<IndexType> iaArray( n + 1, ia );
        const hmemo::HArrayRef<IndexType> jaArray( numValues, ja );
        const hmemo::HArrayRef<ValueType> valueArray( numValues, values );
        setCSRData( rowDist, colDist, numValues, iaArray, jaArray, valueArray );
    }

    /** This method sets raw DIA data in the same way as setDIAData but with raw value array */

    template<typename ValueType>
    void setRawDIAData(
        dmemo::DistributionPtr rowDist,
        dmemo::DistributionPtr colDist,
        const IndexType numDiagonals,
        const IndexType* offsets,
        const ValueType* values )
    {
        //const IndexType numRows    = rowDist->getLocalSize(); // not needed
        const IndexType numColumns = colDist->getLocalSize();
        // use of HArrayRef instead of HArray avoids additional copying of values
        const hmemo::HArrayRef<IndexType> offsetArray( numDiagonals, offsets );
        const hmemo::HArrayRef<ValueType> valueArray( numColumns * numDiagonals, values );
        setDIAData( rowDist, colDist, numDiagonals, offsetArray, valueArray );
    }

    /** Setting raw dense data for a replicated matrix, only for convenience. */

    template<typename ValueType>
    void setRawDenseData( const IndexType n, const IndexType m, const ValueType* values, const ValueType eps = 0 )
    {
        setRawDenseData( dmemo::DistributionPtr( new dmemo::NoDistribution( n ) ),
                         dmemo::DistributionPtr( new dmemo::NoDistribution( m ) ),
                         values,
                         eps );
    }

    /** @brief Assignment of a matrix to this matrix
     *
     * Assignment of a matrix to this matrix with automatic conversion
     * to the matrix type and the current row / column distribution of *this
     *
     *  Note: other.getNumRows() == getNumRows(), other.getNumColumns() == getNumColumns() is mandatory
     *
     *  @param[in] other   the matrix to be converted.
     */
    virtual void assign( const Matrix& other ) = 0;

    /** @brief Assignment of a transposed matrix to this matrix
     *
     *  @param[in] other   the matrix to be assigned.
     */
    virtual void assignTranspose( const Matrix& other ) = 0;

    /**
     * @brief Setting (distributed) matrix with any replicated/global matrix storage data.
     *
     * @param[in] other   replicated matrix data containing all values to be set.
     *
     * Note: Size of other matrix must be exactly the same as this matrix. This routine might imply type and storage format
     * conversions as well as distributing the data according to the current distribution of this matrix.
     * The current matrix data will be overridden.
     */
    virtual void assign( const _MatrixStorage& other ) = 0;

    /** @brief Assignment of local storage that fits a given row distribution.
     *
     *  The columns of the local storage will be splitted according to the column distribution.
     *
     *  @param[in] storage   local part of the matrix on this processor
     *  @param[in] rowDist   the given row distribution.
     *  @param[in] colDist   storage will be splitted according to the column distribution.
     */
    virtual void assign( const _MatrixStorage& storage, dmemo::DistributionPtr rowDist, dmemo::DistributionPtr colDist ) = 0;

    /** @brief Gets the local part (no splitted columns) of a matrix as if the distribution of columns is replicated.
     *
     *  @param[out] storage  will contain the local part of the matrix with all columns.
     *
     *  As splitting of columns is differently handled for sparse and dense matrices,
     *  this method is useful for conversion between them. An alternative solution
     *  of copying the whole matrix and replication of columns might be too expensive.
     */
    virtual void buildLocalStorage( _MatrixStorage& storage ) const = 0;

    /** @brief Get access to the local storage of this matrix according to row and column distribution.
     *
     *  Derived classes might return covariant types.
     */

    virtual const _MatrixStorage& getLocalStorage() const = 0;

    /** @brief This method allows any arbitrary redistribution of the matrix.
     *
     *  @param[in] rowDistribution is new distribution of rows, global size must be mNumRows
     *  @param[in] colDistribution is new distribution of columns, global size must be mNumColumns
     */
    virtual void redistribute( dmemo::DistributionPtr rowDistribution, dmemo::DistributionPtr colDistribution ) = 0;

    /** @brief This method returns one row of the matrix.
     *
     * @param[out] row              is a replicated vector with all values of the row
     * @param[in]  globalRowIndex   global index of the row that should be extracted
     *
     * - the vector row might be of any type but for efficiency it should have the same type as the matrix
     *   (otherwise conversion)
     * - the output vector will always be replicated
     */
    virtual void getRow( Vector& row, const IndexType globalRowIndex ) const = 0;

    /** @brief This method returns the diagonal.
     *
     * @param[out]   diagonal is the destination array
     *
     * Calculations are dependent to the diagonal property.
     */
    virtual void getDiagonal( Vector& diagonal ) const = 0;

    /** @brief This method replaces the diagonal.
     *
     * @param[in] diagonal  is the source array
     *
     * Calculations are dependent to the diagonal property.
     */
    virtual void setDiagonal( const Vector& diagonal ) = 0;

    /** @brief This method replaces the diagonal by a diagonal value.
     *
     * @param[in] scalar  is the source value
     *
     * Calculations are dependent to the diagonal property.
     */
    virtual void setDiagonal( const Scalar scalar ) = 0;

    /** @brief This method scales all values with a vector.
     *
     * @param[in] scaling   is the source array
     *
     * row wise calculations.
     */
    virtual void scale( const Vector& scaling ) = 0;

    /** @brief This method scales all matrix values with a scalar.
     *
     * @param[in] scaling   is the source value.
     */
    virtual void scale( const Scalar scaling ) = 0;

    /** @brief This method replaces all elements with its conjugate value. */

    virtual void conj() = 0;

    /**
     * @brief Returns a copy of the value at the passed global indexes.
     *
     * @param[in] i   the global row index
     * @param[in] j   the global column index
     * @return        a copy of the value at the passed global position.
     *
     * As this operator requires communication in SPMD mode it can be very inefficient in some situations.
     */
    Scalar operator()( IndexType i, IndexType j ) const;

    /**
     * @brief Returns a copy of the value at the passed global indexes.
     *
     * @param[in] i   the global row index
     * @param[in] j   the global column index
     * @return        a copy of the value at the passed global position.
     *
     * As this operation requires communication in SPMD mode it can be very inefficient in some situations.
     */
    virtual Scalar getValue( IndexType i, IndexType j ) const = 0;

    virtual void writeAt( std::ostream& stream ) const;

    /**
     * @brief Returns the number of global rows.
     *
     * @return the number of rows.
     */
    inline IndexType getNumRows() const;

    /**
     * @brief Returns the number of columns.
     *
     * @return the number of columns.
     */
    inline IndexType getNumColumns() const;

    /** @brief Gets the total number of non-zero values in the matrix.
     *
     *  An element is considered to be non-zero if its absolute value
     *  is greater equal than mEpsilon. Zero diagonal elements are also
     *  counted if this->hasDiagonalProperty() is given.
     *
     *  This routine does not count zero elements even if they are stored
     *  (e.g. for dense or dia storage data).
     *
     * @return the number of non-zero values in this matrix.
     */
    virtual IndexType getNumValues() const = 0;

    /**
     * @brief Returns the sparsity rate of the matrix.
     *
     * @return the sparsity rate ( numValues / (numRows * numColumns)) of the whole matrix.
     */
    double getSparsityRate() const;

    /**
     * @brief Returns whether the matrix is symmetric or not.
     *
     * @return a boolean pointing out whether the matrix is symmetric or not.
     */
    bool checkSymmetry() const;

    /**
     * @brief Computes result = alpha * this * x + beta * y.
     *
     * @param[out]  result  the Vector to store the result to
     * @param[in]   alpha   the Scalar alpha of the expression
     * @param[in]   x       the Vector x of the expression
     * @param[in]   beta    the Scalar beta of the expression
     * @param[in]   y       the Vector y of the expression
     *
     * This method computes result = alpha * this * x + beta * y. If
     * result == x or result == y new storage is allocated to store the result.
     */
    virtual void matrixTimesVector(
        Vector& result,
        const Scalar alpha,
        const Vector& x,
        const Scalar beta,
        const Vector& y ) const = 0;

    virtual void vectorTimesMatrix(
        Vector& result,
        const Scalar alpha,
        const Vector& x,
        const Scalar beta,
        const Vector& y ) const = 0;

    /**
     * @brief Computes this = alpha * other.
     *
     * @param[out]  other   the Matrix to multiply
     * @param[in]   alpha   the Scalar of the expression
     */
    virtual void matrixTimesScalar( const Matrix& other, const Scalar alpha ) = 0;

    /**
     * @brief Computes this = alpha * A + beta * B.
     *
     * @param[in]   alpha   the Scalar alpha of the expression
     * @param[in]   A       the Matrix A of the expression
     * @param[in]   beta    the Scalar beta of the expression
     * @param[in]   B       the Matrix B of the expression
     */
    virtual void matrixPlusMatrix( const Scalar alpha, const Matrix& A, const Scalar beta, const Matrix& B ) = 0;

    /**
     * @brief Computes result = alpha * this * B + beta * C.
     *
     * @param[out]  result  the Matrix to store the result to
     * @param[in]   alpha   the Scalar alpha of the expression
     * @param[in]   B       the Matrix B of the expression
     * @param[in]   beta    the Scalar beta of the expression
     * @param[in]   C       the Matrix C of the expression
     */
    virtual void matrixTimesMatrix(
        Matrix& result,
        const Scalar alpha,
        const Matrix& B,
        const Scalar beta,
        const Matrix& C ) const = 0;

    /** Getter routine for the local number of stored values. */

    virtual IndexType getLocalNumValues() const = 0;

    /** Getter routine for the local number of rows. */

    virtual IndexType getLocalNumRows() const = 0;

    /** Getter routine for the local number of columns. */

    virtual IndexType getLocalNumColumns() const = 0;

    /**
     * @brief Gets a constant reference to the column distribution.
     *
     * @return a constant reference to the column distribution.
     */
    inline const dmemo::Distribution& getColDistribution() const;

    /**
     * @brief Gets a pointer to the column distribution.
     *
     * @return a pointer to the column distribution.
     */
    inline dmemo::DistributionPtr getColDistributionPtr() const;

    /**
     * @brief Gets a constant reference to the row distribution.
     *
     * @return a constant reference to the row distribution.
     *
     * Note: same as getDistribution()
     */
    inline const dmemo::Distribution& getRowDistribution() const;

    /**
     * @brief Gets a pointer to the row distribution.
     *
     * @return a pointer to the row distribution.
     *
     * Note: same as getDistribution()
     */
    inline dmemo::DistributionPtr getRowDistributionPtr() const;

    /**
     * @brief Specifies on which compute back end the matrix operations should take place.
     *
     * @param[in] context  the compute back to use for calculations with matrix
     *
     * Note: Only for sparse matrices it is possible to specify separate locations for
     *       local and halo computations.
     */
    virtual void setContextPtr( const hmemo::ContextPtr context ) = 0;

    /**
     * @brief Set individual context for local and halo part of the matrix.
     *
     * @param[in] localContext   context for local part
     * @param[in] haloContext    context for non-local part
     *
     *  Note: Only sparse matrices will override this method, others will ignore second argument.
     */
    virtual void setContextPtr( const hmemo::ContextPtr localContext, const hmemo::ContextPtr haloContext );

    /**
     *  @brief Getter routine for the context.
     *
     *  Note: Only for SparseMatrix the context of the halo can be queried.
     */
    virtual hmemo::ContextPtr getContextPtr() const = 0;

    /**
     * @brief SyncKind describes if the communication and computation should be done synchronously or asynchronously.
     */
    typedef enum
    {
        ASYNCHRONOUS, // asynchronous execution to overlap computations, communications
        SYNCHRONOUS // synchronous, operations will not overlap
    } SyncKind;

    /**
     * @brief MatrixKind describes if a matrix is dense or sparse.
     */
    typedef enum
    {
        DENSE, //!< matrix kind for a dense matrix
        SPARSE //!< matrix kind for a sparse matrix
    } MatrixKind;

    /** Each derived matrix must give info about its kind (DENSE or SPARSE). */

    virtual MatrixKind getMatrixKind() const = 0;

    /**
     * @brief Getter routine for the communication kind.
     *
     * @return the communication kind.
     */
    inline SyncKind getCommunicationKind() const;

    /**
     * @brief Setter routine for the communication kind.
     *
     * @param[in] communicationKind the communication kind.
     */
    void setCommunicationKind( SyncKind communicationKind );

    /**
     * @brief Inherits context and kind arguments from another matrix.
     *
     * @param[in] other   is the input matrix.
     *
     * This routine will also be used by copy constructors in base classes.
     *
     */
    void inheritAttributes( const Matrix& other );

    /**
     * @brief Prefetches matrix data to its 'preferred' context location.
     */
    virtual void prefetch() const = 0;

    /**
     * @brief Waits for a possibly running prefetch.
     */
    virtual void wait() const = 0;

    /**
     * @brief The assignment operator for matrix.
     *
     * @param[in] other   is the input matrix.
     *
     * The assignment operator will make a deep copy of the input matrix. Size and distributions
     * are inherited, but there might be implicit conversions regarding storage format and/or
     * value type of the matrix elements.
     */
    Matrix& operator=( const Matrix& other );

    /**
     * @brief Assignment operator for alhpa * A
     *
     * @param[in] exp   representation of alpha * A as Expression object
     */
    Matrix& operator=( const Expression_SM& exp );

    /**
     * @brief Assignment operator for alhpa * A * B with A and B matrices and scalar alpha
     *
     * @param[in] exp   representation of alpha * A * B as Expression object
     */
    Matrix& operator=( const Expression_SMM& exp );

    /**
     * @brief The assignment operator for a GEMM expression alpha * A * B + beta * C
     *
     * @param[in] exp   representation of alpha * A * B + beta * C as Expression object
     */
    Matrix& operator=( const Expression_SMM_SM& exp );

    /**
     * @brief The assignment operator for alpha * A + beta * B
     *
     * @param[in] exp   expression of the form alpha * A + beta * B
     */
    Matrix& operator=( const Expression_SM_SM& exp );

    /**
     * @brief The assignment operator this *= alpha
     *
     * @param[in] val   Factor used for scaling of the matrix
     */
    Matrix& operator*=( const Scalar val );

    /**
     * @brief The assignment operator this += A
     *
     * @param[in] exp   Matrix to be added
     */
    Matrix& operator+=( const Matrix& exp );

    /**
     * @brief The assignment operator this += alpha * A
     *
     * @param[in] exp   representation of alpha * A as Expression object
     */
    Matrix& operator+=( const Expression_SM& exp );

    /**
     * @brief The assignment operator this -= A
     *
     * @param[in] exp   Matrix to be added
     */
    Matrix& operator-=( const Matrix& exp );

    /**
     * @brief The assignment operator this -= alpha * A
     *
     * @param[in] exp   representation of alpha * A as Expression object
     */
    Matrix& operator-=( const Expression_SM& exp );

    /**
     * @brief Computes the inverse of a matrix.
     *
     * @param[in] other   another matrix with the same shape as this matrix
     *
     * Set this matrix with the inverse of another square matrix.
     * This matrix will have afterwards the same distribution as the
     * other matrix.
     */
    virtual void invert( const Matrix& other ) = 0;

    /**
     * @brief Returns the L1 norm of this.
     *
     * @return the L1 norm of this.
     *
     * l1Norm computes the sum of the absolute values of this.
     */
    virtual Scalar l1Norm() const = 0;

    /**
     * @brief Returns the L2 norm of this.
     *
     * @return the L2 norm of this.
     *
     * l2Norm computes the sum of the absolute values of this.
     */
    virtual Scalar l2Norm() const = 0;

    /**
     * @brief Getter routine of the max norm of this matrix.
     *
     * @return the maximal absolute value for elements of this matrix
     */
    virtual Scalar maxNorm() const = 0;

    /**
     * @brief Returns the max norm of ( this - other ).
     *
     * @param[in] other another matrix with the same shape as this matrix
     * @return the max norm of ( this - other )
     *
     * The maximal value is given by the largest difference between two elements
     * at the same position of the matrices.
     *
     * Note: This method is the most general implementation and should be
     *       implemented more efficiently in derived classes.
     */
    virtual Scalar maxDiffNorm( const Matrix& other ) const;

    /**
     * @brief Constructor function which creates a 'zero' matrix of same type as a given matrix.
     *
     * \code
     * void sub( ..., const Matrix& matrix, ...)
     * {
     *     ...
     *     // Create a copy of the input matrix
     *
     *     common::unique_ptr<Matrix> newmatrix ( matrix.newMatrix() );
     *     *newmatrix = matrix;
     *
     *     // Create a unity matrix of same type and same row distribution as matrix
     *
     *     common::unique_ptr<Matrix> newmatrix ( matrix.newMatrix() );
     *     newmatrix->allocate( matrix.getRowDistributionPtr(), matrix.getRowDistributionPtr() );
     *     newmatrix->setIdentity();
     *     ...
     * }
     * \endcode
     *
     * This method is a workaround to call the constructor of a derived matrix class
     * where the derived class is not known at compile time.
     */
    virtual Matrix* newMatrix() const = 0;

    /*
     *  @brief Create a dense vector with same type/context as matrix and same row distribution
     *
     *  Note: this method is for a more convenient use
     */

    Vector* newDenseVector() const
    {
        Vector* v = Vector::getDenseVector( getValueType(), getRowDistributionPtr() );
        v->setContextPtr( getContextPtr() );
        return v;
    }

    /**
     * @brief Constructor function which creates a copy of this matrix.
     *
     * \code
     * common::unique_ptr<Matrix> newmatrix = matrix.copy();
     * // More convenient to use, but exactly same as follows:
     * common::unique_ptr<Matrix> newmatrix = Matrix::create( matrix.getCreateValue() ); *newmatrix = matrix;
     * \endcode
     *
     * This method is a workaround to call the copy constructor of a derived matrix class
     * where the derived class is not known at compile time.
     */
    virtual Matrix* copy() const = 0;

    /**
     * @brief Copy this matrix with a new distribution
     *
     * @param[in] rowDistribution   new distribution of the rows
     * @param[in] colDistribution   new distribution of the columns
     */
    virtual Matrix* copy( dmemo::DistributionPtr rowDistribution, dmemo::DistributionPtr colDistribution ) const;

    /**
     * @brief Queries the keytype to create matrix from factory with same valuetype and storagetype
     */
    MatrixCreateKeyType getCreateValue() const;

    /**
     * @brief Queries the value type of the matrix elements, e.g. DOUBLE or FLOAT.
     */
    virtual common::scalar::ScalarType getValueType() const = 0;

    /**
     * @brief Query the storage format of the used matrix format, e.g. CSR or ELL.
     */
    virtual MatrixStorageFormat getFormat() const = 0;

    /**
     * @brief Query the size of one matrix element
     */
    virtual size_t getValueTypeSize() const = 0;

    /** Returns the diagonalProperty of the local storage.
     *
     * @return if the diagonal property is full filled.
     */
    virtual bool hasDiagonalProperty() const = 0;

    /**
     * @brief Rechecks the storages for their diagonal property.
     *
     * Usually each matrix has a flag that indicates if the diagonal property is given.
     * This makes the query hasDiagonalProperty very efficient. Therefore matrices
     * keep track of this flag for all their operations and reset it if necessary.
     *
     * This routine is only available to have a workaround if matrix or storage data
     * has been modified by the user outside of the class (NOT RECOMMENDED).
     */

    virtual void resetDiagonalProperty() = 0;

    /**
     * @brief Returns the global memory that is allocated to hold this matrix.
     *
     * getMemoryUsage returns the global memory that is allocated to hold this matrix. For a distributed matrix
     * all partitions are summed together.
     *
     * @return the memory consumption of this matrix.
     */
    virtual size_t getMemoryUsage() const = 0;

protected:

    /**
     * @brief Constructs a matrix by corresponding distributions.
     *
     * @param[in] rowDistribution   specifies number of rows and how the rows of the matrix are
     *                              distributed among the available partitions.
     * @param[in] colDistribution   specifies number of columns and the distribution of a vector
     *                              when matrix is applied to a matrix-vector multiplication.
     *
     *  The colDstribution is used to setup a communication pattern (e.g. halo)
     *  that can be reused for multiple matrix-vector operations. It might also be
     *  used to split the matrix data on one partition into a local and non-local
     *  part.
     *
     *  The number of rows and columns is given implicitly by the global sizes of the distributions.
     */
    Matrix( dmemo::DistributionPtr rowDistribution, dmemo::DistributionPtr colDistribution );

    /**
     * @brief Constructs a matrix of a given size replicated on each partition.
     *
     * @param[in] numRows      number of rows, must be non-negative.
     * @param[in] numColumns   number of columns, must be non-negative.
     *
     * Same as Matrix( NoDistribution(numRows), NoDistribution(numColumns) )
     */
    Matrix( const IndexType numRows, const IndexType numColumns );

    /**
     * @brief Constructs a square matrix with the given size and the specified distribution.
     *
     * @param[in] distribution      specifies how the rows of the matrix are distributed among
     *                              the available partitions.
     *
     * Same as Matrix( distribution, distribution ); column distribution will be same as that of rows.
     */
    Matrix( dmemo::DistributionPtr distribution );

    /**
     * @brief Constructs a square matrix with a replicated distribution.
     *
     * @param[in] size              is the number of rows
     *
     * Same as Matrix(size, size ).
     */
    explicit Matrix( const IndexType size );

    /**
     * @brief Default constructor, creates a replicated matrix of size 0 x 0.
     */
    Matrix();

    /**
     * @brief Redistribute the given matrix with the specified distributions.
     *
     * @param[in] other            the matrix to redistribute
     * @param[in] rowDistribution  specifies how the rows of the matrix are distributed among
     *                             the available partitions.
     * @param[in] colDistribution  specifies the distribution of a vector when matrix is
     *                             applied to a matrix-vector multiplication.
     *
     * The following relations must hold:
     *
     *   - other.getNumRows() == distribution->getGlobalSize()
     *   - other.getNumColumns() == colDistribution->getGlobalSize()
     */
    Matrix( const Matrix& other, dmemo::DistributionPtr rowDistribution, dmemo::DistributionPtr colDistribution );

    /**
     * @brief Copies a matrix to a new matrix with the same distribution.
     *
     * @param[in] other  the matrix to take a copy from
     *
     */
    Matrix( const Matrix& other );

    /**
     * @brief Sets the global/local size of replicated matrix.
     *
     * @param[in] numRows      number of rows, must be non-negative.
     * @param[in] numColumns   number of columns, must be non-negative.
     */
    void setReplicatedMatrix( const IndexType numRows, const IndexType numColumns );

    /**
     * @brief Sets the global and local size of distributed matrix.
     *
     * @param[in] distribution is distribution for rows
     * @param[in] colDistribution is distribution for columns
     *
     * Global and local size is given implicitly by the distributions itself.
     */
    void setDistributedMatrix( dmemo::DistributionPtr distribution, dmemo::DistributionPtr colDistribution );

    void readFromSingleFile( const std::string& fileName );
 
    void readFromPartitionedFile( const std::string& fileName );

    void resetRowDistribution( dmemo::DistributionPtr distribution );

    void resetRowDistributionByFirstColumn();

    dmemo::DistributionPtr mColDistribution;

    // TODO remove mNumRows and mNumColumns, this value is stored in the distribution
    IndexType mNumRows;
    IndexType mNumColumns;

protected:

    void checkSettings() const; // check valid member variables

    void swapMatrix( Matrix& other ); // swap member variables of Matrix

    SCAI_LOG_DECL_STATIC_LOGGER( logger )

private    :

    using Distributed::getDistribution;

    using Distributed::getDistributionPtr;

    void sanityCheck( const Expression<Matrix, Matrix, Times>& exp );

    void sanityCheck( const Expression<Matrix, Matrix, Times>& exp, const Matrix& C );

    void sanityCheck( const Matrix& A, const Matrix& B );

    void setDefaultKind(); // set default values for communication and compute kind

    SyncKind mCommunicationKind;//!< synchronous/asynchronous communication

    void writeToSingleFile(
        const std::string& fileName,
        const std::string& fileType,
        const common::scalar::ScalarType dataType,
        const common::scalar::ScalarType indexType,
        const FileIO::FileMode fileMode ) const;

    void writeToPartitionedFile(
        const std::string& fileName,
        const std::string& fileType,
        const common::scalar::ScalarType dataType,
        const common::scalar::ScalarType indexType,
        const FileIO::FileMode fileMode ) const;
};

/* ======================================================================== */
/*             Inline methods                                               */
/* ======================================================================== */

inline IndexType Matrix::getNumRows() const
{
    //return getRowDistributionPtr().get()->getGlobalSize();
    return mNumRows;
}

inline IndexType Matrix::getNumColumns() const
{
    //return getColDistributionPtr().get()->getGlobalSize();
    return mNumColumns;
}

inline Matrix::SyncKind Matrix::getCommunicationKind() const
{
    return mCommunicationKind;
}

inline const dmemo::Distribution& Matrix::getColDistribution() const
{
    SCAI_ASSERT_ERROR( mColDistribution, "NULL column distribution for Matrix" )
    return *mColDistribution;
}

inline dmemo::DistributionPtr Matrix::getColDistributionPtr() const
{
    SCAI_ASSERT_ERROR( mColDistribution, "NULL column distribution for Matrix" )
    return mColDistribution;
}

inline const dmemo::Distribution& Matrix::getRowDistribution() const
{
    return getDistribution();
}

inline dmemo::DistributionPtr Matrix::getRowDistributionPtr() const
{
    return getDistributionPtr();
}

/** This function prints a SyncKind on an output stream.
 *
 *  \param stream   is the reference to the output stream
 *  \param kind      is the enum value that is printed
 */
inline std::ostream& operator<<( std::ostream& stream, const scai::lama::Matrix::SyncKind& kind )
{
    switch ( kind )
    {
        case scai::lama::Matrix::SYNCHRONOUS:
        {
            stream << "SYNCHRONOUS";
            break;
        }

        case scai::lama::Matrix::ASYNCHRONOUS:
        {
            stream << "ASYNCHRONOUS";
            break;
        }

        default:
        {
            stream << "<unknown sync kind>";
            break;
        }
    }

    return stream;
}

/** This function prints a MatrixKind on an output stream.
 *
 *  \param stream   is the reference to the output stream
 *  \param kind      is the enum value that is printed
 */
inline std::ostream& operator<<( std::ostream& stream, const Matrix::MatrixKind& kind )
{
    switch ( kind )
    {
        case scai::lama::Matrix::DENSE:
        {
            stream << "DENSE";
            break;
        }

        case scai::lama::Matrix::SPARSE:
        {
            stream << "SPARSE";
            break;
        }

        default:
        {
            stream << "<unknown matrix kind>";
            break;
        }
    }

    return stream;
}

} /* end namespace lama */

} /* end namespace scai */
