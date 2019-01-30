/**
 * @file _Matrix.hpp
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
#include <scai/lama/DenseVector.hpp>
#include <scai/lama/io/FileIO.hpp>
#include <scai/lama/matrix/MatrixKind.hpp>
#include <scai/lama/matrix/SyncKind.hpp>

#include <scai/lama/expression/Expression.hpp>
#include <scai/lama/storage/_MatrixStorage.hpp>

// internal scai libraries

#include <scai/common/BinaryOp.hpp>
#include <scai/dmemo/Distribution.hpp>
#include <scai/dmemo/NoDistribution.hpp>

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

typedef std::shared_ptr<class _Matrix> _MatrixPtr;


/** Key type used for the _Matrix factory.
 *
 *  Note: own struct instead of std::pair to allow definition of operator <<
 */

typedef MatrixStorageCreateKeyType MatrixCreateKeyType;

/**
 * @brief The class _Matrix is a abstract type that represents a distributed 2D real or complex matrix.
 *
 * _Matrix is one of the LAMA Base Types and should be used in all situations where it is not necessary to access a
 * single element or to create a new _Matrix.
 */
class COMMON_DLL_IMPORTEXPORT _Matrix:

    public common::Factory<MatrixCreateKeyType, _Matrix*>,
    public dmemo::Distributed
{

public:

    /** @brief More convenient use of the create routine of factory that avoids use of CreateKeyType.
     */
    static _Matrix* getMatrix( const Format format, const common::ScalarType valueType );

    /**
     * @brief Destructor, releases all allocated resources.
     */
    virtual ~_Matrix();

    /**
     * @brief write the matrix to an output file
     *
     * @param[in] fileName is the name of the output file (suffix must be added according to the file type)
     * @param[in] fileType format of the output file (SAMG, _MatrixMarket), default is to decide by suffix
     * @param[in] dataType representation type for output values, default is same type as matrix values
     * @param[in] indexType representation type for col/row index values
     * @param[in] fileMode can be used to forche BINARY or FORMATTED output
     */

    void writeToFile(
        const std::string& fileName,
        const std::string& fileType = "",
        const common::ScalarType dataType = common::ScalarType::UNKNOWN,
        const common::ScalarType indexType = common::ScalarType::UNKNOWN,
        const FileMode fileMode = FileMode::DEFAULT  ) const;

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
     *  @brief Virtual method that delivers a unique identification of the class name to which a matrix belongs.
     *
     *  This method is mainly used for error messages but might also be used to verify
     *  correctness when a base class pointer/reference is statically casted to a derived class.
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
     * void sub( ..., _Matrix& a, ... )
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
     *  \endcode
     */
    void setIdentity( const IndexType n );

    /**
     * This method sets a matrix by reading its values from one or multiple files.
     *
     * @param[in] fileName      the filename to read from
     * @param[in] rowDist       distribution of the matrix
     *
     *   \code
     *      CSRSparseMatrix<double> matrix;
     *      matrix.readFromFile( "matrix.mtx", rowDist )           ! each processor gets its local part of the matrix in one file
     *      matrix.readFromFile( "matrix_%r.mtx", rowDist )        ! read a partitioned matrix with the given distribution
     *   \endcode
     */
    void readFromFile( const std::string& fileName, dmemo::DistributionPtr rowDist );

    /**
     * This method sets a matrix by reading its values from one or multiple files.
     *
     * @param[in] fileName      the filename to read from
     * @param[in] comm          specifies the involved processes
     *
     *   \code
     *      CSRSparseMatrix<double> matrix;
     *      matrix.readFromFile( "matrix.mtx" )       ! read matrix only on processor 0, gets a SingleDistribution
     *      matrix.readFromFile( "matrix_%r.mtx" )    ! general block distributed matrix, each processor reads it own file
     *   \endcode
     */
    void readFromFile( const std::string& fileName, dmemo::CommunicatorPtr comm = dmemo::Communicator::getCommunicatorPtr() );

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

    /** @brief Assignment of a matrix to this matrix
     *
     * Assignment of a matrix to this matrix with automatic conversion
     * to the matrix type and the current row / column distribution of *this
     *
     *  Note: other.getNumRows() == getNumRows(), other.getNumColumns() == getNumColumns() is mandatory
     *
     *  @param[in] other   the matrix to be converted.
     */
    virtual void assign( const _Matrix& other ) = 0;

    /** @brief Assignment of a transposed matrix to this matrix
     *
     *  @param[in] other   the matrix to be assigned.
     */
    virtual void assignTranspose( const _Matrix& other ) = 0;

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
    virtual void assignDistribute( const _MatrixStorage& storage, dmemo::DistributionPtr rowDist, dmemo::DistributionPtr colDist ) = 0;

    /** @brief assignmnet of a matrix with redistribution 
     *
     *  The global size of the new row and column distribution must fit the existing sizes of the input matrix.
     *
     *  @param[in] other     the matrix that will be assigned and redistributed
     *  @param[in] rowDist   the new row distribution
     *  @param[in] colDist   the new column distribution
     */
    virtual void assignDistribute( const _Matrix& other, dmemo::DistributionPtr rowDist, dmemo::DistributionPtr colDist ) = 0;

    virtual void assignLocal( const _MatrixStorage& storage, dmemo::DistributionPtr rowDist ) = 0;

    /** @brief Gets the local part (no splitted columns) of a matrix as if the distribution of columns is replicated.
     *
     *  @param[out] storage  will contain the local part of the matrix with all columns.
     *
     *  The output storage might be of any value type and any format as this operation supports also
     *  implicit type conversion.
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
     *  @param[in] rowDistributionPtr is new distribution of rows, global size must be getNumRows()
     *  @param[in] colDistributionPtr is new distribution of columns, global size must be getNumColumns()
     *
     *  For sparse matrices it might be allowed that the new number of columns might become larger with the
     *  new column distribution. 
     */
    virtual void redistribute( dmemo::DistributionPtr rowDistributionPtr, dmemo::DistributionPtr colDistributionPtr ) = 0;

    /**
     *  @brief Redistribute this matrix with a redistributor 
     *
     *  @param[in] redistributor      object that supports redistribution of rows from old to new distribution
     *  @param[in] colDistributionPtr is new distribution of columns, global size must be getNumColumns()
     *
     *  Note: redistributor.getSourceDistribution() == this->getRowDistribution() must be valid
     *        this->getDistribution() == redistributor.getTargetDistribution() is valid after the call.
     *
     *  This call is more efficient than redistribute( redistributor.getTargetDistributionPtr(), colDistribution )
     *  as the communication schedule is already available and can be reused.
     */
    virtual void redistribute( const dmemo::RedistributePlan& redistributor, dmemo::DistributionPtr colDistributionPtr ) = 0;

    /**
     *  @brief Redistribute this matrix with a redistribution, columns are eiter replicated or have same dist as rows
     */
    virtual void redistribute( const dmemo::RedistributePlan& redistributor );

    /**
     *  @brief Resize a given matrix.
     *
     *  This method is very similiar to redistribute but the new distributions might have a different global size.
     *  Therefore the matrix might be either truncated or filled up with zero in each of its dimension.
     */
    virtual void resize( dmemo::DistributionPtr rowDist, dmemo::DistributionPtr colDist ) = 0;

    /** @brief This method replaces all elements with its conjugate value. */

    virtual void conj() = 0;

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
     *  @brief Getter routine for the context.
     */
    virtual hmemo::ContextPtr getContextPtr() const = 0;

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
    void inheritAttributes( const _Matrix& other );

    /**
     * @brief Prefetches matrix data to its 'preferred' context location.
     */
    virtual void prefetch() const = 0;

    /**
     * @brief Waits for a possibly running prefetch.
     */
    virtual void wait() const = 0;

    /**
     * @brief Computes the inverse of a matrix.
     *
     * @param[in] other   another matrix with the same shape as this matrix
     *
     * Set this matrix with the inverse of another square matrix.
     * This matrix will have afterwards the same distribution as the
     * other matrix.
     */
    virtual void invert( const _Matrix& other ) = 0;

    /**
     * @brief Constructor function which creates a 'zero' matrix of same type as a given matrix.
     *
     * \code
     * void sub( ..., const _Matrix& matrix, ...)
     * {
     *     ...
     *     // Create a copy of the input matrix
     *
     *     std::unique_ptr<Matrix> newmatrix ( matrix.newMatrix() );
     *     *newmatrix = matrix;
     *
     *     // Create a unity matrix of same type and same row distribution as matrix
     *
     *     std::unique_ptr<Matrix> newmatrix ( matrix.newMatrix() );
     *     newmatrix->allocate( matrix.getRowDistributionPtr(), matrix.getRowDistributionPtr() );
     *     newmatrix->setIdentity();
     *     ...
     * }
     * \endcode
     *
     * This method is a workaround to call the constructor of a derived matrix class
     * where the derived class is not known at compile time.
     */
    virtual _Matrix* newMatrix( void ) const = 0;

    /**
     * @brief Constructor function which creates a copy of this matrix.
     *
     * \code
     * _MatrixPtr newmatrix = matrix.copy();
     * // More convenient to use, but exactly same as follows:
     * _MatrixPtr newmatrix = _Matrix::create( matrix.getCreateValue() ); *newmatrix = matrix;
     * \endcode
     *
     * This method is a workaround to call the copy constructor of a derived matrix class
     * where the derived class is not known at compile time.
     */
    virtual _Matrix* copy() const = 0;

    /**
     * @brief Copy this matrix with a new distribution
     *
     * @param[in] rowDistribution   new distribution of the rows
     * @param[in] colDistribution   new distribution of the columns
     */
    virtual _Matrix* copyRedistributed( dmemo::DistributionPtr rowDistribution, dmemo::DistributionPtr colDistribution ) const;

    /**
     * @brief Queries the keytype to create matrix from factory with same valuetype and storagetype
     */
    MatrixCreateKeyType getCreateValue() const;

    /**
     * @brief Queries the value type of the matrix elements, e.g. DOUBLE or FLOAT.
     */
    virtual common::ScalarType getValueType() const = 0;

    /**
     * @brief Query the storage format of the used matrix format, e.g. CSR or ELL.
     */
    virtual Format getFormat() const = 0;

    /**
     * @brief Query the size of one matrix element
     */
    virtual size_t getValueTypeSize() const = 0;

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
    _Matrix( dmemo::DistributionPtr rowDistribution, dmemo::DistributionPtr colDistribution );

    /**
     * @brief Constructs a matrix of a given size replicated on each partition.
     *
     * @param[in] numRows      number of rows, must be non-negative.
     * @param[in] numColumns   number of columns, must be non-negative.
     *
     * Same as _Matrix( NoDistribution(numRows), NoDistribution(numColumns) )
     */
    _Matrix( const IndexType numRows, const IndexType numColumns );

    /**
     * @brief Constructs a square matrix with the given size and the specified distribution.
     *
     * @param[in] distribution      specifies how the rows of the matrix are distributed among
     *                              the available partitions.
     *
     * Same as _Matrix( distribution, distribution ); column distribution will be same as that of rows.
     */
    _Matrix( dmemo::DistributionPtr distribution );

    /**
     * @brief Constructs a square matrix with a replicated distribution.
     *
     * @param[in] size              is the number of rows
     *
     * Same as _Matrix(size, size ).
     */
    explicit _Matrix( const IndexType size );

    /**
     * @brief Default constructor, creates a replicated matrix of size 0 x 0.
     */
    _Matrix();

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
    _Matrix( const _Matrix& other, dmemo::DistributionPtr rowDistribution, dmemo::DistributionPtr colDistribution );

    /**
     * @brief Copies a matrix to a new matrix with the same distribution.
     *
     * @param[in] other  the matrix to take a copy from
     *
     */
    _Matrix( const _Matrix& other );

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

    /** 
     *   @brief Assignment operator for base class _Matrix is explicitly disabled.
     * 
     *   \code
     *       void copyIt( _Matrix& a, const _Matrix& b )
     *       {
     *           a = b;           // NOT allowed
     *           a.assign( b );   // can be used as workaround
     *       }
     *
     *       template<typename ValueType>
     *       void copyIt( Matrix<ValueType>& a, const Matrix<ValueType>& b )
     *       {
     *           a = b;           // This is supported
     *       }
     *   \endcode
     *
     *   Note: assignment operations with expressions are also not supported for the base class _Matrix.
     */
    _Matrix& operator=( const _Matrix& other ) = delete;

    /**
     *   @brief Move assignment operator for base class _Matrix is explicitly disabled.
     *
     *   \code
     *       _Matrix& a = ...
     *       _Matrix& b = ...
     *       a = std::move( b );           // NOT allowed
     *   \endcode
     *
     *   Note: move assignment operator is supported for full matrix classes.
     */
    _Matrix& operator=( _Matrix&& other ) = delete;

    void readFromSingleFile( const std::string& fileName, 
                             dmemo::CommunicatorPtr comm = dmemo::Communicator::getCommunicatorPtr()  );

    void readFromSingleFile( const std::string& fileName, dmemo::DistributionPtr distribution );

    void readFromPartitionedFile( const std::string& fileName, 
                                  dmemo::CommunicatorPtr comm = dmemo::Communicator::getCommunicatorPtr() );

    void resetRowDistribution( dmemo::DistributionPtr distribution );

    void checkSettings() const; // check valid member variables

    void swapMatrix( _Matrix& other ); // swap member variables of _Matrix

    SCAI_LOG_DECL_STATIC_LOGGER( logger )

protected:

    /**
     *  @brief Provide move assignment operator for base class, disable default one 
     *
     *  This method moves only member variables of this base class.
     */
    void moveImpl( _Matrix&& other );

    /**
     *  @brief Check that size of local storage fits the distributions
     */
    static void checkLocalStorageSizes( const  _MatrixStorage& localStorage, const dmemo::Distribution& rowDist );

private:

    /* ============================================================= */
    /*    new member variables of class _Matrix                      */
    /* ============================================================= */

    dmemo::DistributionPtr mColDistribution;

    SyncKind mCommunicationKind;//!< synchronous/asynchronous communication

    using Distributed::getDistribution;

    using Distributed::getDistributionPtr;

    void setDefaultKind(); // set default values for communication and compute kind

    static SyncKind getDefaultSyncKind();  // get default kind as set by environment

    void writeToSingleFile(
        const std::string& fileName,
        const std::string& fileType,
        const common::ScalarType dataType,
        const common::ScalarType indexType,
        const FileMode fileMode ) const;

    void writeToPartitionedFile(
        const std::string& fileName,
        const std::string& fileType,
        const common::ScalarType dataType,
        const common::ScalarType indexType,
        const FileMode fileMode ) const;
};

/* ======================================================================== */
/*             Inline methods                                               */
/* ======================================================================== */

inline IndexType _Matrix::getNumRows() const
{
    return getDistribution().getGlobalSize();
}

inline IndexType _Matrix::getNumColumns() const
{
    return mColDistribution->getGlobalSize();
}

inline SyncKind _Matrix::getCommunicationKind() const
{
    return mCommunicationKind;
}

inline const dmemo::Distribution& _Matrix::getColDistribution() const
{
    return *mColDistribution;
}

inline dmemo::DistributionPtr _Matrix::getColDistributionPtr() const
{
    return mColDistribution;
}

inline const dmemo::Distribution& _Matrix::getRowDistribution() const
{
    return getDistribution();
}

inline dmemo::DistributionPtr _Matrix::getRowDistributionPtr() const
{
    return getDistributionPtr();
}

} /* end namespace lama */

} /* end namespace scai */
