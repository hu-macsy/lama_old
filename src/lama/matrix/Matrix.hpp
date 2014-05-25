/**
 * @file Matrix.hpp
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
 * @brief Abstract base class for all matrices supported by LAMA.
 * @author Jiri Kraus, Thomas Brandes
 * @date 22.02.2011
 * @since 1.0.0
 */
#ifndef LAMA_MATRIX_HPP_
#define LAMA_MATRIX_HPP_

// for dll_import
#include <lama/config.hpp>

// base classes
#include <lama/Distributed.hpp>

// others
#include <lama/LAMATypes.hpp>
#include <lama/Scalar.hpp>
#include <lama/Vector.hpp>
#include <lama/Context.hpp>

#include <lama/distribution/Distribution.hpp>
#include <lama/distribution/NoDistribution.hpp>
#include <lama/expression/Expression.hpp>

// logging
#include <logging/logging.hpp>

namespace lama
{

// Forward declaration

class _MatrixStorage;

/** Pointer class for a matrix, always use of a shared pointer. */

typedef boost::shared_ptr<class Matrix> MatrixPtr;

/**
 * @brief The class Matrix is a abstract type that represents a distributed 2D real or complex matrix.
 *
 * Matrix is one of the LAMA Base Types and should be used in all situations where it is not necessary to access a
 * single element or to create a new Matrix.
 */
class LAMA_DLL_IMPORTEXPORT Matrix: public Distributed
{

public:

    /**
     * @brief ExpressionMemberType is the type that is used the template Expression to store a Vector.
     */
    typedef const Matrix& ExpressionMemberType;

    /**
     * @brief Destructor, releases all allocated resources.
     */
    virtual ~Matrix();

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
     * LAMA_ASSERT_DEBUG( a.isConsistent(), a << ": is invalid matrix after reading" )
     * \endcode
     */

    virtual bool isConsistent() const = 0;

    /**
     *  @brief Virtual method that delivers the class name to which a matrix belongs.
     */
    virtual const char* getTypeName() const = 0;

    /** @brief Creates a LAMA array with the same value type as the matrix.
     *
     *  @return an auto pointer to the LAMA array.
     *
     *  Same as _LAMAArray::create( this.getValueType() )
     *
     *  Value type is known only at runtime, so pointer to the base class
     *  is returned. Auto pointer indicates that calling routine takes ownership of
     *  the allocated array.
     */
    _LAMAArray* createArray() const;

    /**
     * @brief Clears the full matrix, resets global and local sizes to 0.
     *
     * \code
     *     CSRSparseMatrix<double> a ( ... )
     *     a = CSRSparseMatrix<double> ();     \\ will free all arrays
     *     a.clear();                          \\ same functionality, clears involved arrays
     *
     * \endcode
     */
    virtual void clear() = 0;

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
    virtual void allocate( DistributionPtr rowDistribution, DistributionPtr colDistribution ) = 0;

    /**
     * @brief Operator that sets the matrix to the identity matrix.
     *
     * \code
     * void sub( ..., Matrix& a, ... )
     * ...
     * LAMA_ASSERT_EQUAL_DEBUG( a.getNumRows(), a.getNumColumns() )
     * a.setIdentity( a.getRowDistribution() );
     * \endcode
     */

    /** Set matrix to a identity square matrix with same row and column distribution. */

    virtual void setIdentity( DistributionPtr distribution ) = 0;

    /** Set matrix to a (replicated) identity matrix with same row and column distribution. */

    void setIdentity( const IndexType n );

    /**
     * This method sets a matrix by reading its values from a file.
     *
     * @param[in] filename      the filename to read from
     *
     * Each matrix class must provide an implementation of this method.
     * The matrix might have any distribution.
     */
    virtual void readFromFile( const std::string& filename ) = 0;

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
    virtual void setDenseData( 
        DistributionPtr rowDist,
        DistributionPtr colDist, 
        const _LAMAArray& values,
        Scalar eps = Scalar( 0 )        ) = 0;

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
        DistributionPtr rowDist, 
        DistributionPtr colDist, 
        const IndexType numValues,
        const LAMAArray<IndexType>& ia, 
        const LAMAArray<IndexType>& ja, 
        const _LAMAArray& values ) = 0;

    /** This method sets raw dense data in the same way as setDenseData but with raw value array */

    template<typename ValueType>
    void setRawDenseData( 
        DistributionPtr rowDist, 
        DistributionPtr colDist, 
        const ValueType* values,
        const ValueType eps = 0 )
    {
        const IndexType n = rowDist->getLocalSize();
        const IndexType m = colDist->getGlobalSize();

        // use of LAMAArrayRef instead of LAMAArray avoids additional copying of values

        const LAMAArrayRef<ValueType> valueArray( values, n * m );

        setDenseData( rowDist, colDist, valueArray, Scalar( eps ) );
    }

    /** This method sets raw CSR data in the same way as setCSRData but with raw value array */

    template<typename ValueType>
    void setRawCSRData(
        DistributionPtr rowDist,
        DistributionPtr colDist,
        const IndexType numValues,
        const IndexType* ia,
        const IndexType* ja,
        const ValueType* values )
    {
        const IndexType n = rowDist->getLocalSize();

        // use of LAMAArrayRef instead of LAMAArray avoids additional copying of values

        const LAMAArrayRef<IndexType> iaArray( ia, n + 1 );
        const LAMAArrayRef<IndexType> jaArray( ja, numValues );
        const LAMAArrayRef<ValueType> valueArray( values, numValues );

        setCSRData( rowDist, colDist, numValues, iaArray, jaArray, valueArray );
    }

    /** Setting raw dense data for a replicated matrix, only for convenience. */

    template<typename ValueType>
    void setRawDenseData( 
        const IndexType n, 
        const IndexType m, 
        const ValueType* values,
        const ValueType eps = 0 )
    {
        setRawDenseData( DistributionPtr( new NoDistribution( n ) ),
                         DistributionPtr( new NoDistribution( m ) ), 
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
    virtual void assign( const _MatrixStorage& storage, DistributionPtr rowDist, DistributionPtr colDist ) = 0;

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
    virtual void redistribute( DistributionPtr rowDistribution, DistributionPtr colDistribution ) = 0;

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
    inline const Distribution& getColDistribution() const;

    /**
     * @brief Gets a pointer to the column distribution.
     *
     * @return a pointer to the column distribution.
     */
    inline DistributionPtr getColDistributionPtr() const;

    /**
     * @brief Specifies on which compute back end the matrix operations should take place.
     *
     * @param[in] context  the compute back to use for calculations with matrix
     *
     * Note: Only for sparse matrices it is possible to specify separate locations for
     *       local and halo computations.
     */
    virtual void setContext( const ContextPtr context ) = 0;

    /**
     * @brief Set individual context for local and halo part of the matrix.
     *
     * @param[in] localContext   context for local part
     * @param[in] haloContext    context for non-local part
     *
     *  Note: Only sparse matrices will override this method, others will ignore second argument.
     */
    virtual void setContext( const ContextPtr localContext, const ContextPtr haloContext );

    /**
     *  @brief Getter routine for the context.
     *
     *  Note: Only for SparseMatrix the context of the halo can be queried.
     */
    virtual ContextPtr getContextPtr() const = 0;

    /**
     * @brief Method returns a reference to the constant context.
     *
     * @return    reference to the constant context.
     */
    virtual const Context& getContext() const
    {
        return *getContextPtr();
    }

    /**
     * @brief SyncKind describes if the communication and computation should be done synchronously or asynchronously.
     */
    typedef enum
    {
        ASYNCHRONOUS,   // asynchronous execution to overlap computations, communications
        SYNCHRONOUS     // synchronous, operations will not overlap
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
     */
    virtual Scalar maxDiffNorm( const Matrix& other ) const = 0;

    /**
     * @brief Constructor function which creates a 'zero' matrix of same type as a given matrix.
     *
     * \code
     * void sub( ..., const Matrix& matrix, ...)
     * {
     *     ...
     *     // Create a copy of the input matrix
     *
     *     std::auto_ptr<Matrix> newMatrix ( matrix.create() );
     *     *newMatrix = matrix;
     *
     *     // Create a unity matrix of same type and same row distribution as matrix
     *
     *     std::auto_ptr<Matrix> newMatrix ( matrix.create() );
     *     newMatrix->allocate( matrix.getRowDistributionPtr(), matrix.getRowDistributionPtr() );
     *     newMatrix->setIdentity();
     *     ...
     * }
     * \endcode
     *
     * This method is a workaround to call the constructor of a derived matrix class
     * where the derived class is not known at compile time.
     */
    virtual Matrix* create() const = 0;

    /**
     * @brief Constructor creates a distributed zero matrix of same type as a given matrix.
     *
     * @param[in] rowDistribution   TODO[doxy] Complete Description.
     * @param[in] colDistribution   TODO[doxy] Complete Description.
     */
    Matrix* create( DistributionPtr rowDistribution, DistributionPtr colDistribution ) const;

    /**
     * @brief Constructor creates a distributed dense vector of same type as a given matrix.
     *
     * @param[in] distribution   TODO[doxy] Complete Description.
     * @param[in] value          TODO[doxy] Complete Description.
     */
    Vector* createDenseVector( DistributionPtr distribution, const Scalar value ) const;

    /**
     * @brief Constructor function which creates a copy of this matrix.
     *
     * \code
     * std::auto_ptr<Matrix> newMatrix = matrix.copy();
     * // More convenient to use, but exactly same as follows:
     * std::auto_ptr<Matrix> newMatrix = matrix.create(); *newMatrix = matrix;
     * \endcode
     *
     * This method is a workaround to call the copy constructor of a derived matrix class
     * where the derived class is not known at compile time.
     */
    virtual Matrix* copy() const = 0;

    /**
     * @brief Queries the value type of the matrix elements, e.g. DOUBLE or FLOAT.
     */
    virtual Scalar::ScalarType getValueType() const = 0;

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
    Matrix( DistributionPtr rowDistribution, DistributionPtr colDistribution );

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
    Matrix( DistributionPtr distribution );

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
    Matrix( const Matrix& other, DistributionPtr rowDistribution, DistributionPtr colDistribution );

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
    void setDistributedMatrix( DistributionPtr distribution, DistributionPtr colDistribution );

    DistributionPtr mColDistribution;

    // TODO remove mNumRows and mNumColumns, this value is stored in the distribution
    IndexType mNumRows;
    IndexType mNumColumns;

protected:

    void checkSettings() const; // check valid member variables

    void swapMatrix( Matrix& other ); // swap member variables of Matrix

    LAMA_LOG_DECL_STATIC_LOGGER( logger )

private:

    void sanityCheck( const Expression<Matrix, Matrix, Times>& exp );

    void sanityCheck( const Expression<Matrix, Matrix, Times>& exp, const Matrix& C );

    void sanityCheck( const Matrix& A, const Matrix& B );

    void setDefaultKind(); // set default values for communication and compute kind

    SyncKind mCommunicationKind;//!< synchronous/asynchronous communication
};

/* ======================================================================== */
/*             Inline methods                                               */
/* ======================================================================== */

inline IndexType Matrix::getNumRows() const
{
    //return getDistributionPtr().get()->getGlobalSize();
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

inline const Distribution& Matrix::getColDistribution() const
{
    LAMA_ASSERT_ERROR( mColDistribution, "NULL column distribution for Matrix" )
    return *mColDistribution;
}

inline DistributionPtr Matrix::getColDistributionPtr() const
{
    LAMA_ASSERT_ERROR( mColDistribution, "NULL column distribution for Matrix" )
    return mColDistribution;
}

/** This function prints a SyncKind on an output stream.
 *
 *  \param stream   is the reference to the output stream
 *  \param kind      is the enum value that is printed
 */
inline std::ostream& operator<<( std::ostream& stream, const Matrix::SyncKind& kind )
{
    switch ( kind )
    {
    case Matrix::SYNCHRONOUS:
    {
        stream << "SYNCHRONOUS";
        break;
    }
    case Matrix::ASYNCHRONOUS:
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
    case Matrix::DENSE:
    {
        stream << "DENSE";
        break;
    }
    case Matrix::SPARSE:
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

}

#endif // LAMA_MATRIX_HPP_
