/**
 * @file Matrix.hpp
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
 * @brief Abstract base class for all matrices of a certain value type.
 * @author Thomas Brandes
 * @date 31.10.2017
 */
#pragma once

#include <scai/lama/matrix/_Matrix.hpp>

#include <scai/lama/matrix/MatrixAssembly.hpp>
#include <scai/lama/freeFunction.hpp>

#include <scai/lama/storage/MatrixStorage.hpp>
#include <scai/lama/expression/CastMatrixExpression.hpp>
#include <scai/lama/expression/ComplexMatrixExpression.hpp>
#include <scai/lama/expression/MatrixExpressions.hpp>

namespace scai
{

namespace lama
{

/**
 * @brief Abstract base class for a matrix of a certain value type.
 */
template<typename ValueType>
class COMMON_DLL_IMPORTEXPORT Matrix:

    public _Matrix

{

public:

    /**
     * @brief ExpressionMemberType is the type that is used to represent a matrix in a template expression.
     *
     * When a matrix is used in an expression, a const reference is used. 
     */
    typedef const Matrix<ValueType>& ExpressionMemberType;

    /**
     * @brief Define ObjectValueType so matrix class can be used in certain free functions to deduce ValueType.
     */
    typedef ValueType ObjectValueType;

    /** Create a new matrix of a certain format but with same value type */

    static Matrix<ValueType>* getMatrix( const Format format );

    virtual ~Matrix();

    /** Overwrite _Matrix::newMatrix to get the covariant return type */

    virtual Matrix<ValueType>* newMatrix() const = 0;

    /** 
     *  This method returns a new allocated target vector ( has getRowDistribution() )
     *
     *  This method might be used in polymorphic code where the result vector of
     *  matrixTimesVector is not always a dense vector.
     */
    virtual Vector<ValueType>* newTargetVector() const;

    /** 
     *  This method returns a new allocated source vector ( has getColDistribution() )
     */
    virtual Vector<ValueType>* newSourceVector() const;

    /** Overwrite _Matrix::copy to get the covariant return type */

    virtual Matrix<ValueType>* copy( void ) const = 0;

    /** Overwrite _Matrix::getLocalStorage to get the covariant return type */

    virtual const MatrixStorage<ValueType>& getLocalStorage( void ) const = 0;

    /**
     * @brief The assignment operator for matrix.
     *
     * @param[in] other   is the input matrix.
     *
     * The assignment operator will make a deep copy of the input matrix. Size and distributions
     * are inherited, but there might be implicit conversions regarding storage format.
     */
    Matrix<ValueType>& operator=( const Matrix<ValueType>& other );

    /**
     * @brief Assign of a matrix with an implicit operation, e.g. transpose( A ), conj( A )
     *
     * @param[in] other   matrix and implicit operation
     *
     * \code
     *     DenseMatrix<double> m( DenseStorage<double>( 2, 3, HArray<double>( { 1, 2, 3, 4, 5, 6} ) ) );
     *     DenseMatrix<double> mT;
     *     mT = transpose( m );
     * \endcode
     */
    Matrix<ValueType>& operator=( const OpMatrix<ValueType>& other );

    /**
     * @brief The assignment operator with explicit type conversion.
     *
     * \code
     *     Matrix<float>& fM = ...
     *     const Matrix<double>& dM = ...
     *     fM = cast<float>( dM );
     * \endcode
     */
    template<typename OtherValueType>
    Matrix<ValueType>& operator=( const CastMatrixExpression<ValueType, OtherValueType>& exp );

    /** this = real( A ) or this = imag( A ) */

    template<common::ComplexPart kind, typename OtherValueType>
    Matrix<ValueType>& operator=( const ComplexPartMatrixExpression<OtherValueType, kind>& exp );

    /** this = cmplx( A, B ) */

    Matrix<ValueType>& operator=( const ComplexBuildMatrixExpression<RealType<ValueType> >& exp );

    /**
     * @brief Assignment operator for alhpa * A
     *
     * @param[in] exp   representation of alpha * A as Expression object
     */
    Matrix<ValueType>& operator=( const Expression_SM<ValueType>& exp );

    /**
     * @brief Assignment operator for alhpa * A * B with A and B matrices and scalar alpha
     *
     * @param[in] exp   representation of alpha * A * B as Expression object
     */
    Matrix<ValueType>& operator=( const Expression_SMM<ValueType>& exp );

    /**
     * @brief The assignment operator for a GEMM expression alpha * A * B + beta * C
     *
     * @param[in] exp   representation of alpha * A * B + beta * C as Expression object
     */
    Matrix<ValueType>& operator=( const Expression_SMM_SM<ValueType>& exp );

    /**
     * @brief The assignment operator for alpha * A + beta * B
     *
     * @param[in] exp   expression of the form alpha * A + beta * B
     */
    Matrix<ValueType>& operator=( const Expression_SM_SM<ValueType>& exp );

    /**
     * @brief The assignment operator this += A
     *
     * @param[in] exp   Matrix to be added, must have same type
     */
    Matrix& operator+=( const Matrix<ValueType>& exp );

    /**
     * @brief The assignment operator this -= A
     *
     * @param[in] exp   Matrix to be added, must have same type
     */
    Matrix& operator-=( const Matrix<ValueType>& exp );

    /**
     * @brief The assignment operator this += alpha * A
     *
     * @param[in] exp   representation of alpha * A as Expression object
     */
    Matrix<ValueType>& operator+=( const Expression_SM<ValueType>& exp );

    /**
     * @brief The assignment operator this -= alpha * A
     *
     * @param[in] exp   representation of alpha * A as Expression object
     */
    Matrix<ValueType>& operator-=( const Expression_SM<ValueType>& exp );

    /**
     *  @brief Use of assignment operator '*=' to scale this matrix.
     */
    Matrix<ValueType>& operator*=( const ValueType alpha );
    /**
     * @brief Returns a copy of the value at the passed global indexes.
     *
     * @param[in] i   the global row index
     * @param[in] j   the global column index
     * @return        a copy of the value at the passed global position.
     *
     * As this operator requires communication in SPMD mode it can be very inefficient in some situations.
     */
    ValueType operator()( IndexType i, IndexType j ) const;

    /**
     * @brief Assign this matrix a diagonal matrix specified by a vector containing the diagonal elements.
     *
     * @param[in] diagonal contains the values for the diagonal
     *
     * \code
     *     Matrix<ValueType>& m = ...
     *     // m.setIdentity( n ) can also be written as follows
     *     m.assignDiagonal( fill<DenseVector<ValueType>>( n, ValueType( 1 ) ) );  
     * \endcode
     */
    virtual void assignDiagonal( const Vector<ValueType>& diagonal ) = 0;

    /**
     * @brief Returns a copy of the value at the passed global indexes.
     *
     * @param[in] i   the global row index
     * @param[in] j   the global column index
     * @return        a copy of the value at the passed global position.
     *
     * As this operation requires communication in SPMD mode it can be very inefficient in some situations.
     */
    virtual ValueType getValue( IndexType i, IndexType j ) const = 0;

    /**
     * @brief Update of an (existing ) element in a matrix
     *
     * @param[in] i   the global row index
     * @param[in] j   the global column index
     * @param[in] val value used for update
     * @param[in] op  binary operation used to combine new and old value, default is COPY
     *
     * Note: this method will never change the pattern of a sparse matrix.
     */
    virtual void setValue(
        const IndexType i,
        const IndexType j,
        const ValueType val,
        const common::BinaryOp op = common::BinaryOp::COPY ) = 0;

    /** 
     *  Implementation of pure method of _Matrix::getValueType 
     *
     *  The following code demonstrates how this routine might be used to
     *  make a safe static_cast from the base class _Matrix.
     *
     *  \code
     *    _Matrix& m = ...   
     *    if ( m.getValueType() == TypeTraits<T>::stype )
     *    {
     *        Matrix<T>& mt = static_cast<Matrix<T>&>( m );
     *    }
     *  \endcode
     */
    virtual common::ScalarType getValueType() const;

    /** Implementation of pure method _Matrix::getValueTypeSize */

    virtual size_t getValueTypeSize() const;

    /**
     * @brief Returns whether the matrix is symmetric or not.
     *
     * @return a boolean pointing out whether the matrix is symmetric or not.
     */
    bool checkSymmetry() const;

    /**
     * @brief Computes result = alpha * op( this ) * x + beta * y.
     *
     * @param[out]  result  the vector to store the result to
     * @param[in]   alpha   the scalar alpha of the expression
     * @param[in]   x       the vector x of the expression
     * @param[in]   beta    the scalar beta of the expression
     * @param[in]   y       the vector y of the expression, optional
     * @param[in]   op      operation (transpose, conj) implicitly applied to this matrix
     *
     */
    virtual void matrixTimesVector(
        Vector<ValueType>& result,
        const ValueType alpha,
        const Vector<ValueType>& x,
        const ValueType beta,
        const Vector<ValueType>* y,
        const common::MatrixOp op ) const;

    /**
     *  @brief Special case of matrixTimesVector but with all vectors are dense 
     *
     * @param[out] result        dense vector that stores the result
     * @param[in]  alpha         scaling factor for matrix * vector
     * @param[in]  x             dense vector that is used for multiplication
     * @param[in]  beta          scaling factor for additional summand
     * @param[in]  y             additional summand ( beta = 0 if not available )
     * @param[in]  op            specifies the matrix operation (transpose, conj) implicitly applied
     *
     * This virtual method must be overridden by all derived classes that have no
     * own version of matrixTimesVector.
     */
    virtual void matrixTimesVectorDense(
        DenseVector<ValueType>& result,
        const ValueType alpha,
        const DenseVector<ValueType>& x,
        const ValueType beta,
        const DenseVector<ValueType>* y,
        const common::MatrixOp op ) const = 0;

    /**
     * @brief Computes this = alpha * A.
     *
     * @param[out]  A       the matrix to multiply
     * @param[in]   alpha   the Scalar of the expression
     */
    virtual void matrixTimesScalar( const Matrix<ValueType>& A, const ValueType alpha ) = 0;

    /** Elementwise binary operation of matrix elements
     *
     *  @param[in] matrixA, matrixB are the input matrices, must have the same distribution
     *  @param[in] op               specifies the binary operation to be applied
     *
     *  This matrix becomes the result of the operation, alias with one of the input matrices is supported.
     */
    virtual void binaryOp( const Matrix<ValueType>& matrixA, const common::BinaryOp op, const Matrix<ValueType>& matrixB ) = 0;

    /**
     * @brief Computes this = alpha * A + beta * B.
     *
     * @param[in]   alpha   the Scalar alpha of the expression
     * @param[in]   A       the _Matrix A of the expression
     * @param[in]   beta    the Scalar beta of the expression
     * @param[in]   B       the _Matrix B of the expression
     */
    virtual void matrixPlusMatrix( const ValueType alpha, const Matrix<ValueType>& A, 
                                   const ValueType beta,  const Matrix<ValueType>& B ) = 0;

    /**
     * @brief Computes result = alpha * this * B + beta * C.
     *
     * @param[out]  result  the matrix to store the result to
     * @param[in]   alpha   the scalar alpha of the expression
     * @param[in]   B       the matrix B of the expression
     * @param[in]   beta    the scalar beta of the expression
     * @param[in]   C       the matrix C of the expression
     */
    virtual void matrixTimesMatrix(
        Matrix<ValueType>& result,
        const ValueType alpha,
        const Matrix<ValueType>& B,
        const ValueType beta,
        const Matrix<ValueType>& C ) const = 0;

    /** Implementation of _Matrix::setRow for all typed matrices
     *
     */
    void setRow( const _Vector& row, const IndexType globalRowIndex,
                 const common::BinaryOp op );

    /** Implementation of _Matrix::setColumn for all typed matrices
     *
     *  The method is implemented by setting the local part of the column on each partition.
     *  All derived classes must provide setLocalColum( colArray, colIndex, op )
     */
    void setColumn( const _Vector& column,
                    const IndexType colIndex,
                    const common::BinaryOp op );

    /**
     * @brief Returns the L1 norm of this matrix.
     *
     * @return the L1 norm of this as real value 
     *
     * l1Norm computes the sum of the absolute values of all entries
     */
    virtual RealType<ValueType> l1Norm( void ) const = 0;

    /**
     * @brief Returns the L2 norm of this.
     *
     * @return the L2 norm of this.
     *
     * l2Norm computes the sum of the absolute values of this.
     */
    virtual RealType<ValueType> l2Norm( void ) const = 0;

    /**
     * @brief Returns the max norm of this matrix
     *
     * @return the maximal absolute value for elements of this matrix
     */
    virtual RealType<ValueType> maxNorm( void ) const = 0;

    /**
     * @brief Returns the max norm of ( this - other ).
     *
     * @param[in] other another matrix with the same shape as this matrix
     * @return the max norm of ( this - other )
     *
     * The maximal value is given by the largest difference between two elements
     * at the same position of the matrices. This method must be implemented by
     * derived classes.
     */
    virtual RealType<ValueType> maxDiffNorm( const Matrix<ValueType>& other ) const = 0;

    /* ======================================================================= */
    /*     setter / getter for diagonal of a matrix                            */
    /* ======================================================================= */

    /** @brief Get the diagonal of a (square) matrix
     *
     * @param[out]   diagonal will contain the diagonal of this matrix
     *
     *  This matrix must be a square matrix with the same row and column distribution.
     *  Note: diagonal will have the same distribution.
     */
    virtual void getDiagonal( Vector<ValueType>& diagonal ) const = 0;

    /** @brief This method replaces the diagonal
     *
     * @param[in] diagonal  contains the new diagonal, must have row distribution of matrix
     *
     * For a sparse matrix, the matrix must have the diagonal property, i.e. the sparse
     * pattern has an entry for each diagonal element.
     */
    virtual void setDiagonal( const Vector<ValueType>& diagonal ) = 0;

    /** @brief This method replaces the diagonal by a diagonal value.
     *
     * @param[in] scalar  is the source value
     *
     * Calculations are dependent to the diagonal property.
     */
    virtual void setDiagonal( const ValueType& scalar ) = 0;

    /* ======================================================================= */
    /*     scaling of matrix entries                                           */
    /* ======================================================================= */

    /** @brief This method scales all matrix elements with a scalar value.
     *
     * @param[in] alpha is the scaling factor.
     */
    virtual void scale( const ValueType& alpha ) = 0;

    /** @brief This method scales the matrix elements individually for each row.
     *
     * @param[in] scaleY  is a vector whose distribution must match the row distribution
     *
     * This operation corresponds to $this = diagonalMatrix( scaleY ) * this$, i.e.
     * pre-multiplying this matrix with a diagonal marix built by the vector scaleY.
     */
    virtual void scaleRows( const DenseVector<ValueType>& scaleY ) = 0;

    /** @brief This method scales the matrix elements individually for each column.
     *
     * @param[in] scaleY  is a vector whose distribution must match the column distribution
     *
     * This operation corresponds to $this = this * diagonalMatrix( scaleY )$, i.e.
     * post-multiplying this matrix with a diagonal marix built by the vector scaleY.
     */
    virtual void scaleColumns( const DenseVector<ValueType>& scaleY ) = 0;

    /* ======================================================================= */
    /*     set/get of rows/columns of a matrix                                 */
    /* ======================================================================= */

    /** @brief This method returns one row of this matrix.
     *
     * @param[out] row              is the vector that will contain the queried row of this matrix
     * @param[in]  globalRowIndex   global index of the row that should be extracted
     *
     * - The result vector will have the same distribution as the column distribution of this matrix
     * - the vector row might be of any value type but for efficiency it should have the same value type as this matrix
     * - row might be a sparse or a dense vector, but it is recommended to use a dense vector to get the row 
     *   of a dense matrix and a sparse vector on a sparse matrix.
     * - This method implies communication as one row resides only on one processor but here the communication
     *   pattern for sparse or dense matrices are exploited.
     */
    virtual void getRow( Vector<ValueType>& row, const IndexType globalRowIndex ) const = 0;

    /** @brief This method returns a row of this matrix locally for one processor.
     *
     * @param[out] row            is the vector that will contain the queried row of this matrix
     * @param[in]  localRowIndex  local index of the row that should be extracted
     *
     * - The result vector is not distributed, only valid results on the corresponding partition
     * - the vector row might be of any value type but for efficiency it should have the same value type as this matrix
     * - row might be a sparse or a dense vector, but it is recommended to use a dense vector to get the row 
     *   of a dense matrix and a sparse vector on a sparse matrix.
     * - This method is completely local, no communication
     *   pattern for sparse or dense matrices are exploited.
     */
    virtual void getRowLocal( Vector<ValueType>& row, const IndexType localRowIndex ) const = 0;

    /** @brief This method returns one column of the matrix.
     *
     * @param[out] column           is a distributed vector with all values of the col
     * @param[in]  globalColIndex   global column index of the col that should be extracted
     *
     * - the vector column might be of any type but for efficiency it should have the same type as the matrix
     *   (otherwise conversion)
     * - the distribution of col will be the same as the row distribution of the matrix
     */
    virtual void getColumn( Vector<ValueType>& column, const IndexType globalColIndex ) const = 0;

    /** @brief This method sets one row of the matrix.
     *
     * @param[in]  row              is a non-distributed vector
     * @param[in]  globalRowIndex   global row index of the row that should be set
     * @param[in]  op               specifies the binary op how to combine old and new element
     *
     * - the vector row might be of any type but for efficiency it should have the same type as the matrix
     *   (otherwise conversion)
     * - this method throws an exception for a sparse matrix if the pattern must be changed
     */
    virtual void setRow( const Vector<ValueType>& row,
                         const IndexType globalRowIndex,
                         const common::BinaryOp op );

    /** @brief Method to set one column of the matrix.
     *
     * @param[in]  column           is a distributed vector with all values of the col
     * @param[in]  globalColIndex   global column index of the col that should be set
     * @param[in]  op               specifies the binary op how to combine old and new element
     *
     * - the distribution of col should be the same as the row distribution of the matrix (otherwise temporary)
     * - this method does not change the pattern of a sparse matrix, so throws an exception if it is insufficient
     *
     *  Note: all derived classes must provide setLocalRow( rowArray, localRowIndex, op )
     */
    virtual void setColumn(
        const Vector<ValueType>& column,
        const IndexType globalColIndex,
        const common::BinaryOp op );

    /* ======================================================================= */
    /*     Reductions for matrices                                             */
    /* ======================================================================= */

    /** @brief This method reduces the rows ( dim = 0 ) to a column vector or the columns ( dim = 1 ) 
     *         to a row vector.
     *
     *  @param[out] v is the result vector, has rowDistribution for dim = 1 or colDistribution for dim = 0 
     *  @param[in]  dim must be either 0 or 1
     *  @param[in]  reduceOp specifies operation used for reduction, e.g. ADD, MIN, MAX
     *  @param[in]  elemOp specfies operatin applied to the elements before reduction
     *
     *  \code
     *     const _Matrix& m; _Vector& v;
     *     m.reduce( v, dim = 0, common::BinaryOp::ADD, common::BinaryOp::SQR );  // builds row sums 
     *     m.reduce( v, dim = 1, common::BinaryOp::ADD, common::BinaryOp::SQR );  // builds diagonal of m' m 
     *  \endcode
     */

    virtual void reduce(
        Vector<ValueType>& v,
        const IndexType dim,
        const common::BinaryOp reduceOp,
        const common::UnaryOp elemOp ) const = 0;

    /* ======================================================================= */
    /*     set/fill assembly matrix                                            */
    /* ======================================================================= */

    /** Fill this matrix with assembled matrix data.
     *
     *  @param[in] assembly contains the assembled entries (individually by each processor)
     *  @param[in] op       either COPY or ADD, specifies how to deal with entries at same positions
     *
     *  The matrix must already have beeen allocated before this method is called.
     */
    virtual void fillFromAssembly( const MatrixAssembly<ValueType>& assembly, common::BinaryOp op = common::BinaryOp::COPY );

    /**
     *  @brief Insert all non-zero coefficients of this matrix in an assembly
     *
     *  Each processor adds its part locally with global coordinates.
     *
     *  @param[in,out] assembly    is the object to which non-zero coefficients of this matrix are inserted
     *  @param[in]     rowOffset   is a global offset that is added to the row coordinates
     *  @param[in]     colOffset   is a global offset that is added to the column coordinates
     *
     *  The offsets can be used to reposition this matrix when it is combined/joined with other matrices.
     */
    virtual void disassemble( 
        MatrixAssembly<ValueType>& assembly, 
        const IndexType rowOffset = 0,
        const IndexType colOffset = 0 ) const = 0;

    /* ======================================================================= */
    /*     Concatenation of matrices                                           */
    /* ======================================================================= */

    /**
     * @brief Concatenate multiple matrices horizontally/vertically to a new matrix.
     *
     * @param[in] rowDist   specifies the distribution of the rows for the concatenated matrix
     * @param[in] colDist   specifies the distribution of the columns for the concatenated matrix
     * @param[in] matrices  variable number of const references/pointers to the matrices
     *
     * The routine decides by its arguments how the matrices will be concatenated. As the size of 
     * the result matrix is explicitly specified, the input matrices are row-wise filled up.
     *
     * This routine should also be able to deal with aliases, i.e. one ore more of the input matrices  might be
     * the pointer to the result matrix.
     */
    virtual void concatenate( 
        dmemo::DistributionPtr rowDist, 
        dmemo::DistributionPtr colDist, 
        const std::vector<const Matrix<ValueType>*>& matrices );

    /**
     *   shorthand for cat( 0, { &m1, &m2 }, 2 )
     */
    void vcat( const Matrix<ValueType>& m1, const Matrix<ValueType>& m2 );

    /**
     *   shorthand for cat( 1, { &m1, &m2 }, 2 )
     */
    void hcat( const Matrix<ValueType>& m1, const Matrix<ValueType>& m2 );

    /**
     *  General routine set a matrix globally with dense data 
     */
    void setRawDenseData( const IndexType n, const IndexType m, const ValueType* values );

    Matrix();

    /**
     * @brief Constructs a matrix of a given size replicated on each partition.
     *
     * @param[in] numRows      number of rows, must be non-negative.
     * @param[in] numColumns   number of columns, must be non-negative.
     *
     * Same as _Matrix( NoDistribution(numRows), NoDistribution(numColumns) )
     */
    Matrix( const IndexType numRows, const IndexType numColumns );

    Matrix( dmemo::DistributionPtr rowDistribution, dmemo::DistributionPtr colDistribution );

    Matrix( const _Matrix& other, dmemo::DistributionPtr rowDistribution, dmemo::DistributionPtr colDistribution );

    Matrix( const _Matrix& other );

    Matrix( const Matrix<ValueType>& other );

    virtual void setLocalRow( const hmemo::HArray<ValueType>& row,
                              const IndexType localRowIndex,
                              const common::BinaryOp op  ) = 0;

    virtual void setLocalColumn( const hmemo::HArray<ValueType>& column,
                                 const IndexType colIndex,
                                 const common::BinaryOp op  ) = 0;

    /** This method is the same for dense/sparse matrices as column distribution is replicated */

    void vectorTimesMatrixRepCols(
        DenseVector<ValueType>& denseResult,
        const ValueType alphaValue,
        const DenseVector<ValueType>& denseX,
        const ValueType betaValue,
        const DenseVector<ValueType>* denseY ) const;

public:

    /** 
     *  @brief This method selects the real or imaginay part of a complex matrix
     *
     *  @param[out] x    is the matrix that will contain the real or imaginary part of this matrix
     *  @param[in]  kind specifies which part (real or complex) is selected
     */
    virtual void selectComplexPart( Matrix<RealType<ValueType> >& x, common::ComplexPart kind ) const = 0;

    /** 
     *  @brief This method builds this complex matrix by two matrices, one contains the real parts, the other the imaginary parts.
     *
     *  @param[in] x is the matrix containing the real parts
     *  @param[in] y is the matrix containing the complex parts
     */
    virtual void buildComplex( const Matrix<RealType<ValueType> >& x, const Matrix<RealType<ValueType> >& y ) = 0;
};

/** 
 * Definiton of corresponding shared pointer type for the class Matrix<ValueType> by a type alias.
 *
 *  \code
 *      MatrixPtr<ValueType> x( Matrix<ValueType>::getMatrix( Format::COO ) );
 *      std::shared_ptr<Matrix<ValueType> > x( Matrix<ValueType>::getMatrix( Format::COO ) );
 *  \endcode
*/
template<typename ValueType>
using MatrixPtr = std::shared_ptr<Matrix<ValueType> >;

/** 
 * Definiton of corresponding unique pointer type for the class Matrix<ValueType> by a type alias.
 *
 *  \code
 *      MatrixPtr1<ValueType> x( Matrix<ValueType>::getMatrix( Format::COO ) );
 *      std::unique_ptr<Matrix<ValueType> > x( Matrix<ValueType>::getMatrix( Format::COO ) );
 *  \endcode
*/
template<typename ValueType>
using MatrixPtr1 = std::unique_ptr<Matrix<ValueType> >;

/* ======================================================================================= */
/*     Implementation of inline methods                                                    */
/* ======================================================================================= */

template<typename ValueType>
template<typename OtherValueType>
Matrix<ValueType>& Matrix<ValueType>::operator=( const CastMatrixExpression<ValueType, OtherValueType>& exp )
{
    this->assign( exp.getArg() );
    return *this;
}

template<typename ValueType>
template<common::ComplexPart kind, typename OtherValueType>
Matrix<ValueType>& Matrix<ValueType>::operator=( const ComplexPartMatrixExpression<OtherValueType, kind>& exp )
{
    // use a static assert to check for correct types of method selectComplexPart, otherwise strange error messages

    static_assert( std::is_same<ValueType, RealType<OtherValueType> >::value,
                   "realMatrix = real|imag( complexMatrix ), value type of realMatrix is not real type of complexMatrix" );

    exp.getArg().selectComplexPart( *this, kind );
    return *this;
}

template<typename ValueType>
Matrix<ValueType>& Matrix<ValueType>::operator=( const ComplexBuildMatrixExpression<RealType<ValueType> >& exp )
{
    buildComplex( exp.getRealArg(), exp.getImagArg() );
    return *this;
}


} /* end namespace lama */

} /* end namespace scai */

