/**
 * @file Vector.hpp
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
 * @brief Definition of an abstract class for distributed vectors of a given type
 * @author Thomas Brandes
 * @date 30.10.2017
 */
#pragma once

#include <scai/lama/_Vector.hpp>

#include <scai/lama/VectorAssembly.hpp>
#include <scai/lama/freeFunction.hpp>

#include <scai/lama/expression/UnaryVectorExpression.hpp>
#include <scai/lama/expression/CastVectorExpression.hpp>
#include <scai/lama/expression/ComplexVectorExpression.hpp>
#include <scai/lama/expression/VectorExpressions.hpp>

#include <scai/common/TypeTraits.hpp>

#include <memory>
#include <type_traits>

namespace scai
{

namespace lama
{

/**
 * @brief Definition of an abstract class that represents a distributed one-dimensional vector 
 *        of a certain value type.
 *
 * @tparam ValueType stands for the type of the entries in the distributed vector.
 *
 * Derived classes might support dense, sparse, and joined vectors.
 */
template <typename ValueType>
class COMMON_DLL_IMPORTEXPORT Vector: public _Vector
{
public:

    /**
     * @brief ExpressionMemberType is the type that is used the template Expression to store a _Vector.
     */
    typedef const Vector<ValueType>& ExpressionMemberType;

    /**
     * @brief Define ObjectValueType so Vector class can be used in certain free functions.
     */
    typedef ValueType ObjectValueType;

    /** Create a new vector of a certain kind but with same value type */

    static Vector<ValueType>* getVector( VectorKind kind );

    /** Desctructor. */

    virtual ~Vector();

    /** Overwrite _Vector::newVector to get the covariant return type */

    virtual Vector<ValueType>* newVector( void ) const = 0;

    /** Overwrite _Vector::copy to get the covariant return type */

    virtual Vector<ValueType>* copy( void ) const = 0;

    /** Implementation of _Vector::getValueType */

    virtual common::ScalarType getValueType() const;

    using _Vector::getContext;
    using _Vector::getDistribution;

    /** Help class to observe the further use of operator[] for Vector */

    class VectorElemProxy
    {
    public:

        /** Proxy constructed by ref to the array and the index value. */

        inline VectorElemProxy( Vector<ValueType>& vector, const IndexType i );

        /** Proxy for a vector element can be used to get its value, type conversion to ValueType
         *
         *  @returns current value of the vector element as a single value 
         */
        inline operator ValueType() const;

        /** indexed value proxy can be assigned a scalar */

        inline VectorElemProxy& operator= ( ValueType val );

        /** Override the default assignment operator to avoid ambiguous interpretation of a[i] = b[i] */

        inline VectorElemProxy& operator= ( const VectorElemProxy& other );

    private:

        Vector<ValueType>& mVector;
        IndexType mIndex;
    };

    /**
     * @brief Returns the value at the passed global index.
     *
     * @param[in] globalIndex   the global index to get the value at.
     * @return                  a copy of the value at the passed global position.
     *
     * As this operation requires communication in SPMD mode it can be very inefficient in some situations.
     * Therefore it is recommended to query values on the local vector data with local indexes.
     */
    virtual ValueType getValue( IndexType globalIndex ) const = 0;

    /**
     * @brief Provide operator() for getValue to make its use more convenient.
     *
     * @param[in] i    the global index to get the value at.
     * @return         a copy of the value at the passed global position.
     *
     * As this operator requires communication ins SPMD mode it can be very inefficient in some situations.
     */
    ValueType operator()( const IndexType i ) const;

    /**
     *
     * @brief This methods sets/updates a value of a vector.
     *
     * Be careful: this method might throw an exception on a sparse vector, if the element is not available
     */
    virtual void setValue( const IndexType globalIndex, const ValueType value ) = 0;

    /**
     *  Indexing of a distributed vector returns a proxy so that this operator can be used
     *  on lhs and rhs of an assignment.
     */
    VectorElemProxy operator[]( const IndexType i );

    /**
     *  Indexing of a const distributed vector returns directly the corresponding element.
     */
    ValueType operator[]( const IndexType i ) const;

    /** 
     * @brief Assignment 'vector *= value' to scale all elements of a vector. 
     * @param[in] value   the value to multiply all elements of this with.
     * @return            a reference to this.
     */
    Vector& operator*=( const ValueType value );

    /** @brief Assignment 'vector /= value' is same as vector.scale( 1 / value ) */

    Vector& operator/=( const ValueType value );

    /**
     * @brief Elementwise multiplication with another vector (same size), i.e. this[i] = this[i] * other[i]
     *
     * @param[in] other   the vector to multiply to do the multiplication per element
     *
     * Note: the other vector can be any type, no temporary is created here
     */
    Vector& operator*=( const Vector<ValueType>& other );

    template<typename OtherValueType>
    Vector& operator*=( const CastVectorExpression<ValueType, OtherValueType>& other );

    /**
     * @brief Elementwise division with another vector (same size), i.e. this[i] = this[i] / other[i]
     *
     * @param[in] other   the vector used for elementwise division
     *
     * Note: the other vector can be any type, no temporary is created here
     */
    Vector& operator/=( const Vector<ValueType>& other );

    template<typename OtherValueType>
    Vector& operator/=( const CastVectorExpression<ValueType, OtherValueType>& other );

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

    /**
     * @brief Returns the max norm of this.
     *
     * @return the max norm of this.
     *
     * maxNorm computes the value of this with the largest magnitude.
     */
    virtual RealType<ValueType> maxNorm() const = 0;

    /**
     * @brief Returns the max norm of the difference with another vector
     *
     *  v1.maxDiffNorm( v2 ) is equivalent to:
     *
     *  \code
     *      Vector<ValueType> tmp = v1 - v2;
     *      maxNorm( tmp )
     *  \endcode
     *
     *  But it avoids the temporary vector wherever possible
     */
    virtual RealType<ValueType> maxDiffNorm( const Vector<ValueType>& other ) const = 0;

    /**
     * @brief Returns the global minimum value of this.
     *
     * @return   the global minimum value of this vector.
     */
    virtual ValueType min() const = 0;

    /**
     * @brief Returns the global maximum value of this.
     *
     * @return the global maximum value of this vector.
     */
    virtual ValueType max() const = 0;

    /**
     * @brief Returns the sum of all vector elements.
     *
     * @return the sum of all vector elements.
     *
     * As the summation of the values depends on the mapping of the values to
     * the processors, this routine might return slightly different results
     * for different parallel environments.
     */
    virtual ValueType sum() const = 0;

    /**
     * @brief Returns the dot product of this and other.
     *
     * @param[in] other   the vector to calculate the dot product with.
     * @return            the dot product of this and other
     */
    virtual ValueType dotProduct( const Vector<ValueType>& other ) const = 0;

    /* =========================================================== */
    /*     this = <vector_expression>                              */
    /* =========================================================== */

    /** this = alpha * A * x */

    Vector<ValueType>& operator=( const Expression_SMV<ValueType>& expression );

    /** this = alpha * x + beta * y */

    Vector<ValueType>& operator=( const Expression_SV_SV<ValueType>& expression );

    /** this = alpha * A * x + beta * y */

    Vector<ValueType>& operator=( const Expression_SMV_SV<ValueType>& expression );

    /** this = alpha * x */

    Vector<ValueType>& operator=( const Expression_SV<ValueType>& expression );

    /** this = alpha * x + beta */

    Vector<ValueType>& operator=( const Expression_SV_S<ValueType>& );

    /** this = x * y */

    Vector<ValueType>& operator=( const Expression_VV<ValueType>& );

    /** this = alpha * x * y */

    Vector<ValueType>& operator=( const Expression_SVV<ValueType>& );

    /** this = x binaryOp y */

    template<common::BinaryOp op>
    Vector<ValueType>& operator=( const Expression<intern::Scalar, Vector<ValueType>, op>& exp );

    template<common::BinaryOp op>
    Vector<ValueType>& operator=( const Expression<Vector<ValueType>, intern::Scalar, op>& exp );

    template<common::BinaryOp op>
    Vector<ValueType>& operator=( const Expression<Vector<ValueType>, Vector<ValueType>, op>& exp );

    /** set same value for all vector elements */

    Vector<ValueType>& operator=( const ValueType value );

    /** @brief Add a scalar elementwise to a vector by using operator $+=$ */

    Vector<ValueType>& operator+=( const ValueType value );

    /** @brief Subtract a scalar elementwise to a vector by using operator $-=$  */

    Vector<ValueType>& operator-=( const ValueType value );

    Vector<ValueType>& operator=( const Vector<ValueType>& other );

    /** @brief this += vector, supported also for mixed value types */

    Vector<ValueType>& operator+=( const Vector<ValueType>& other );

    template<typename OtherValueType>
    Vector<ValueType>& operator+=( const CastVectorExpression<ValueType, OtherValueType>& exp );

    /** @brief this -= vector, vector must have same value type */

    Vector<ValueType>& operator-=( const Vector<ValueType>& other );

    /** @brief this -= cast<ValueType>( vector )  */

    template<typename OtherValueType>
    Vector<ValueType>& operator-=( const CastVectorExpression<ValueType, OtherValueType>& exp );

    /** this +=  alpha * A * x */

    Vector<ValueType>& operator+=( const Expression_SMV<ValueType>& expression );

    /** this +=  alpha * x */

    Vector<ValueType>& operator+=( const Expression_SV<ValueType>& expression );

    /** this -=  alpha * A * x */

    Vector<ValueType>& operator-=( const Expression_SMV<ValueType>& expression );

    /** this -=  alpha * x */

    Vector<ValueType>& operator-=( const Expression_SV<ValueType>& expression );

    /** this = unaryop( x ) */

    Vector<ValueType>& operator=( const UnaryVectorExpression<ValueType>& unaryVectorExp );

    template<typename OtherValueType>
    Vector<ValueType>& operator=( const CastVectorExpression<ValueType, OtherValueType>& exp );

    /** this = real( x ) or this = imag( x ) */

    template<common::ComplexPart kind, typename OtherValueType>
    Vector<ValueType>& operator=( const ComplexPartVectorExpression<OtherValueType, kind>& exp );

    /** this = cmplx( x, y ) */

    Vector<ValueType>& operator=( const ComplexBuildVectorExpression<RealType<ValueType> >& exp );

    /** Initialization of an allocated vector with one value, does not change size/distribution  */
  
    virtual void setScalar( const ValueType& alpha ) = 0;

    /**
     *  @brief Set this vector by applying a unary operation to another vector.
     * 
     *  @param[in] x is the input vector
     *  @param[in] op specifies the unary operation applied to x
     * 
     *  The call y.unaryOp( x, op ) is equivalent to the following code:
     *
     *  \code
     *      for ( IndexType i = 0; i < v1.size(); ++i )
     *      {
     *          y[i] = op( x[i] ); 
     *      }
     *  \endcode
     *
     *  Here are some examples how this method is used:
     *  \code
     *      y = sin( x )  : y.unaryOp( x, common::UnaryOp::SIN )
     *      y = 1 / x     : y.unaryOp( x, common::UnaryOp::RECIPROCAL )
     *  \endcode
     *
     *  Alias of x with this vector is allowed and must be supported, should be faster
     *  due to in-place updates.
     */
    virtual void unaryOp( const Vector<ValueType>& x, const common::UnaryOp op ) = 0;

    /**
     *  @brief Convenient abbreviation for unaryOp where input vector this vector.
     *  
     *  \code
     *     lama::DenseVector<ValueType> x;
     *     ....
     *     x = sin ( x );       // text book syntax
     *     x.unaryOp( x, common::UnaryOp::SIN );
     *     x.unaryOpInPlace( common::UnaryOp::SIN );
     *  \endcode;
     */
    void unaryOpInPlace( const common::UnaryOp op );

    /** Apply binary operation elementwise to elements of vector */

    virtual void binaryOp( const Vector<ValueType>& x, const common::BinaryOp op, const Vector<ValueType>& y ) = 0;

    /**
     *  @brief Set this vector by applying a binary operation to another vector and a scalar
     * 
     *  @param[in] x is the input vector
     *  @param[in] alpha is the scalar element used for the operation
     *  @param[in] op specifies the binary operation for the update
     *  @param[in] swap if true the operands are swapped (op is not always commutative)
     * 
     *  The call y.binaryOp( alhpa, x, op, swap ) is equivalent to the following code:
     *
     *  \code
     *      for ( IndexType i = 0; i < v1.size(); ++i )
     *      {
     *          x[i] = y[i] op s;    // swap = false
     *          x[i] = s op y[i];    // swap = true
     *      }
     *  \endcode
     *
     *  Here are some examples how this method is used:
     *  \code
     *     v = pow( alpha, v )  -> v.binaryOpScalar( v, alpha, common::BinaryOp::POW, true )
     *     v = pow( v, alpha )  -> v.binaryOpScalar( v, alpha, common::BinaryOp::POW, false )
     *     v = 2 / x            -> v.binaryOpScalar( x, ValueType( 2 ), common::BinaryOp::DIVIDE, true )
     *  \endcode
     */
    virtual void binaryOpScalar( const Vector<ValueType>& x, const ValueType& alpha, const common::BinaryOp op, const bool swap ) = 0;

    /** @brief Method for binarOp with scalar as left operand (for convenience) 
     *
     *  \code
     *     DenseVector<ValueType> x( .. );
     *     DenseVector<ValueType> y;
     *     ValueType alpha = 3;
     *     y = alpha / x;     
     *     y.binaryOpScalar( x, alpha, common::BinaryOp::DIVIDE, true );   // internal use only
     *     y.binaryOp( alpha, common::BinaryOp::DIVIDE, x );
     *  \endcode
     */
    inline void binaryOp( const ValueType& alpha, const common::BinaryOp op, const Vector<ValueType>& x );

    /** 
     *  @brief This method selects the real or imaginay part of a complex vector.
     *
     *  @param[out] x    is the vector that will contain either the real or imaginary part of this vector
     *  @param[in]  part specifies whether real or imaginary part is selected
     */
    virtual void selectComplexPart( Vector<RealType<ValueType> >& x, common::ComplexPart part ) const = 0;

    /** 
     *  @brief This method builds a complex vector by two vectors, one contains the real parts ,the other the imaginary parts.
     * 
     *  @param[in] x vector that is taken for the real parts
     *  @param[in] y vector that is taken for the imaginary parts
     *
     *  Like for other binary operations of vectors the input vectors must have same distribution. This
     *  result vector inherits the distributioon from the input vectors.
     *
     *  Note: this operations is not yet supported if one vector is just a scalar, i.e. the real or imaginary 
     *        part is the same for all elements. As workaround you might use a sparse vector with the corresponding 
     *        zero element and no non-zero values.
     */
    virtual void buildComplex( const Vector<RealType<ValueType> >& x, const Vector<RealType<ValueType> >& y ) = 0;

    /** @brief Method for binarOp with scalar as right operand (for convenience) 
     *
     *  \code
     *     DenseVector<ValueType> x( .. );
     *     DenseVector<ValueType> y;
     *     ValueType alpha = 3;
     *     y = pow( x, alpha );                                          // text book syntax
     *     y.binaryOpScalar( x, alpha, common::BinaryOp::POW, false );   // internal use only
     *     y.binaryOp( x, common::BinaryOp::POW, alpha );                // more conveninent method call
     *  \endcode
     */
    inline void binaryOp(  const Vector<ValueType>& x, const common::BinaryOp op, const ValueType& alpha );

    /**
     * @brief Assignment of a 'full' vector expression this = alpha * x + beta * y 
     *
     * Each vector class has to implement its own version of this assignment. 
     */
    virtual void vectorPlusVector( const ValueType& alpha, const Vector<ValueType>& x, 
                                   const ValueType& beta, const Vector<ValueType>& y ) = 0;

    /**
     * @brief Element-wise multiplication of two vectors and a scalar
     *
     * Each vector class has to implement its own version of this assignment. 
     */
    virtual void vectorTimesVector( const ValueType& alpha, const Vector<ValueType>& x, const Vector<ValueType>& y ) = 0;

    /**
     * @brief Assignment of a 'full' vector expression vectorResult = scalarAlpha * vectorX * scalarBeta
     *
     * Each vector class has to implement its own version of this assignment. 
     */
    virtual void vectorPlusScalar( const ValueType& alpha, const Vector<ValueType>& x, const ValueType& beta ) = 0;

    /**
     *  @brief Boolean reduction returns true if all elements fullfill the compare operation with a scalar.
     */
    virtual bool all( common::CompareOp op, const ValueType alpha ) const = 0;

    /**
     *  @brief Boolean reduction returns true if elementwise comparison with other vector is true for all elements
     */
    virtual bool all( common::CompareOp op, const Vector<ValueType>& x ) const = 0;

    /** 
     *  This method gives the vector a size and initializes it with a value.
     *
     *  @param[in] n is the size of the replicated vector
     *  @param[in] value is the value assigned to all elements
     *
     *  \code
     *    DenseVector<double> v1; 
     *    v1.setSameValue( n, value );
     *    DenseVector<double> v2( n );
     *    v2 = value;
     *    DenseVector<double> v3( n, value );
     *  \endcode
     */
    void setSameValue( const IndexType n, const ValueType value )
    {
        allocate( n );
        setScalar( value );
    }

    /** 
     *  This method gives the vector a distribution and initializes it with a value.
     *
     *  @param[in] dist specifies size of the vector and mapping to the processors
     *  @param[in] value is the value assigned to all elements
     *
     *  \code
     *    DenseVector<double> v1; 
     *    v1.setSameValue( n, value );
     *    DenseVector<double> v2( n );
     *    v2 = value;
     *    DenseVector<double> v3( n, value );
     *  \endcode
     */
    void setSameValue( dmemo::DistributionPtr dist, const ValueType value )
    {
        allocate( dist );
        setScalar( value );
    }

    /** Set a replicated vector with sparse vector data
     *
     *  @param[in] n will be the size of the vector
     *  @param[in] nonZeroIndexes positions with non-zero values
     *  @param[in] nonZeroValues values for the non-zero value
     *  @param[in] zeroValue is the 'zero' value, defaults to 0
     *
     *  nonZeroIndexes and nonZeroValues must have the same size. nonZeroIndexes must 
     *  contain valid indexes. They do not have to be sorted.
     */
    void setSparseData(
        const IndexType n,
        const hmemo::HArray<IndexType>& nonZeroIndexes,
        const hmemo::_HArray& nonZeroValues,
        const ValueType zeroValue = ValueType( 0 ) )
    {
        setSameValue( n, zeroValue );
        fillSparseData( nonZeroIndexes, nonZeroValues, common::BinaryOp::COPY );
    }

    /** Fill this vector with assembled vector data.
     *
     *  @param[in] assembly contains the assembled entries (individually by each processor)
     *  @param[in] op       either COPY or ADD, specifies how to deal with entries at same positions
     *
     *  The vector must already have beeen allocated before this method is called.
     */
    void fillFromAssembly( const VectorAssembly<ValueType>& assembly, common::BinaryOp op = common::BinaryOp::COPY );

    /**
     *  @brief Disassemble all vector entries in an assembly
     *
     *  @param[in,out] assembly  is the object to which non-zero values of this vector are added
     *  @param[in]     offset    is a global offset that is added to the coordinates
     *
     *  \code
     *      const Vector<ValueType>& v1 = ...
     *      const Vector<ValueType>& v2 = ...
     *      // concatenate and distribute
     *      VectorAssembly assembly;
     *      v1.disassemble( assembly );
     *      v2.disassemble( assembly, v1.size() );
     *      dist = std::make_shared<BlockDistribution>( v1.size() + v2.size() );
     *      auto v = fill<DenseVector<ValueType>>( dist, 0 );
     *      v.fillFromAssembly( assembly );
     *  \endcode
     */
    void disassemble( VectorAssembly<ValueType>& assembly, const IndexType offset = 0 ) const;

    /** Same as setSparseData but here with raw data for non-zero indexes and values. 
     *
     *  @tparam OtherValueType is the type of the raw data 
     *  @param[in] n will be the size of the vector
     *  @param[in] nnz stands for the number of the non-zero values
     *  @param[in] nonZeroIndexes pointer to array with positions of non-zero values
     *  @param[in] nonZeroValues pointer to array with values
     *  @param[in] zeroValue is the value for all positions that do not appear in nonZeroIndexes
     *
     *  Note: The value type of the raw data might be different to the value type of the vector.
     */
    template<typename OtherValueType>
    void setSparseRawData(
        const IndexType n,
        const IndexType nnz,
        const IndexType nonZeroIndexes[],
        const OtherValueType nonZeroValues[],
        const ValueType zeroValue = ValueType( 0 ) );

    /**
     * This method initilaizes all values of an allocated vector with random numbers.
     *
     * @param[in] bound draw random numbers in the range between 0 and bound (inclusive)
     *
     * For complex vectors a random value is drawn for each real and imaginary part.
     *
     * Keep in mind that bound is an integer value. If you need randonm numbers with other numerical
     * boundaries you should scale them as follows:
     *
     * \code
     *     DistributionPtr dist ( ... );
     *     DenseVector<ValueType> v( dist );
     *     ValueType lb = -1.5, ub = 2.6;
     *     v.fillRandom( 1 );
     *     A = lb + v * ( ub - lb );   // random numbers in the range of lb .. ub
     * \endcode
     */
    virtual void fillRandom( const IndexType bound ) = 0;

    /**
     * @brief Allocate a replicated vector and fill it with random numbers
     *
     * \code
     *    v.setRandom( n, bound ); // same as:  v.allocate( n ); v.fillRandom( bound );
     * \endcode
     *
     * Be careful: in a parallel environment each processor might initialize the array with different
     * values. By calling Math::srandom( seed ) with the same seed on each processor, it can be forced
     * to have the same values.
     */
    void setRandom( const IndexType n, const IndexType bound );

    /**
     * This method sets a distributed vector by its distribution and initializes it with random numbers.
     */
    void setRandom( dmemo::DistributionPtr dist, const IndexType bound );

    /**
     *  Allocate a vector by its size, initialize it with a zero value and fill it sparsely.
     *
     * \code
     *     A.allocate( n );
     *     A = zeroValue;
     *     A.fillSparseRandom( fill, bound );
     * \endcode
     */
    void setSparseRandom( const IndexType n, const ValueType& zeroValue, const float fillRate, const IndexType bound );

    /**
     *  Allocate a vector by its distribution, initialize it with a zero value and fill it sparsely.
     *
     * \code
     *     A.allocate( dist );
     *     A = zeroValue;
     *     A.fillSparseRandom( fill, bound );
     * \endcode
     */
    void setSparseRandom( dmemo::DistributionPtr dist, const ValueType& zeroValue, const float fillRate, const IndexType bound );

    /**
     * @brief Concatenate multiple vectors to a new vector.
     *
     * @param[in] dist specifies the distribution of the concatenated vector.
     * @param[in] vectors is a vector with const pointers/references to the concatenated vectors
     *
     * Note: dist.getGlobalSize() == v[0]->size() + ... v[n-1]->size() 
     *
     * This routine should also be able to deal with aliases, i.e. one ore more of the pointers might be
     * this vector itself.
     */
    virtual void concatenate( dmemo::DistributionPtr dist, const std::vector<const Vector<ValueType>*>& vectors );

    /**
     * @brief Concatenate two vectors to a new vector.
     *
     * @param[in] v1 first part of the new vector
     * @param[in] v2 second part of the new vector
     */
    virtual void cat( const Vector<ValueType>& v1, const Vector<ValueType>& v2 );

protected:

    /**
     * Constructor of replicated vector for derived classes by size and/or context
     *
     * @param[in] size    number of entries for the vector
     * @param[in] context is optional, will be Host context.
     *
     * Note:ithis constructor overrides also the default constructor.
     */
    explicit Vector( const IndexType size = 0, hmemo::ContextPtr context = hmemo::ContextPtr() );

    /**
     * @brief Constructor of Vector for derived classes by distribution
     *
     * @param[in] distribution  the distribution to use for the new Vector.
     * @param[in] context       is optional, will be Host context.
     */
    explicit Vector( dmemo::DistributionPtr distribution, hmemo::ContextPtr context = hmemo::ContextPtr() );

    /**
     * @brief Creates a copy of the passed _Vector.
     *
     * @param[in] other   the Vector to take a copy from.
     *
     * Inherits size/distribution, context and content of the passed vector
     */
    explicit Vector( const _Vector& other );

    /** Override the default copy constructor */

    Vector( const Vector<ValueType>& other );
};
  
/** 
 * Definiton of corresponding shared pointer type for the class Vector<ValueType> by a type alias.
 *
 *  \code
 *      VectorPtr<ValueType> x( Vector<ValueType>::getVector( VectorKind::SPARSE ) );
 *      std::shared_ptr<Vector<ValueType> > x( Vector<ValueType>::getVector( VectorKind::DENSE ) );
 *  \endcode
*/
template<typename ValueType>
using VectorPtr = std::shared_ptr<Vector<ValueType> >;

/** 
 * Definiton of corresponding unique pointer type for the class Vector<ValueType> by a type alias.
 *
 *  \code
 *      VectorPtr1<ValueType> x( Vector<ValueType>::getVector( VectorKind::SPARSE ) );
 *      std::unique_ptr<Vector<ValueType> > x( Vector<ValueType>::getVector( VectorKind::DENSE ) );
 *  \endcode
*/
template<typename ValueType>
using VectorPtr1 = std::unique_ptr<Vector<ValueType> >;

/* ------------------------------------------------------------------------- */
/*  Implementation of basic inline methods                                   */
/* ------------------------------------------------------------------------- */

template<typename ValueType>
ValueType Vector<ValueType>::operator()( const IndexType i ) const
{
    return getValue( i );
}

template<typename ValueType>
ValueType Vector<ValueType>::operator[]( const IndexType i ) const
{
    return getValue( i );
}

template<typename ValueType>
typename Vector<ValueType>::VectorElemProxy Vector<ValueType>::operator[]( const IndexType i )
{
    return VectorElemProxy( *this, i );
}

/* ------------------------------------------------------------------------- */
/*  Implementation of inline methods for VectorElemProxy                     */
/* ------------------------------------------------------------------------- */

template<typename ValueType>
Vector<ValueType>::VectorElemProxy::VectorElemProxy( Vector<ValueType>& vector, const IndexType i ) :

    mVector( vector ),
    mIndex( i )

{
}

template<typename ValueType>
Vector<ValueType>::VectorElemProxy::operator ValueType() const
{
    return mVector.getValue( mIndex );
}

template<typename ValueType>
typename Vector<ValueType>::VectorElemProxy& Vector<ValueType>::VectorElemProxy::operator= ( ValueType val )
{
    mVector.setValue( mIndex, val );
    return *this;
}

template<typename ValueType>
typename Vector<ValueType>::VectorElemProxy& Vector<ValueType>::VectorElemProxy::operator= ( const Vector<ValueType>::VectorElemProxy& other )
{
    ValueType tmp = other.mVector.getValue( other.mIndex );
    mVector.setValue( mIndex, tmp );
    return *this;
}

/* ------------------------------------------------------------------------- */
/*  Implementation of inline methods                                         */
/* ------------------------------------------------------------------------- */

template<typename ValueType>
template<typename OtherValueType>
void Vector<ValueType>::setSparseRawData(
    const IndexType n,
    const IndexType nnz,
    const IndexType nonZeroIndexes[],
    const OtherValueType nonZeroValues[],
    const ValueType zeroValue )
{
    setSameValue( n, zeroValue );
    hmemo::HArrayRef<IndexType> aNonZeroIndexes( nnz, nonZeroIndexes );
    hmemo::HArrayRef<OtherValueType> aNonZeroValues( nnz, nonZeroValues );
    fillSparseData( aNonZeroIndexes, aNonZeroValues, common::BinaryOp::COPY );
}

/* ------------------------------------------------------------------------- */
/*  Inline translation for operators                                         */
/* ------------------------------------------------------------------------- */

template<typename ValueType>
Vector<ValueType>& Vector<ValueType>::operator=( const Vector<ValueType>& other )
{
    setVector( other, common::BinaryOp::COPY );
    return *this;
}

template<typename ValueType>
template<typename OtherValueType>
Vector<ValueType>& Vector<ValueType>::operator=( const CastVectorExpression<ValueType, OtherValueType>& other )
{
    setVector( other.getArg(), common::BinaryOp::COPY );
    return *this;
}

template<typename ValueType>
template<common::ComplexPart kind, typename OtherValueType>
Vector<ValueType>& Vector<ValueType>::operator=( const ComplexPartVectorExpression<OtherValueType, kind>& exp )
{
    // use a static assert to check for correct types of method selectComplexPart, otherwise strange error messages

    static_assert( std::is_same<ValueType, RealType<OtherValueType> >::value, 
                   "realVector = real/imag( complexVector ), value type of realVector is not real type of complexVector" );

    exp.getArg().selectComplexPart( *this, kind );
    return *this;
}

template<typename ValueType>
Vector<ValueType>& Vector<ValueType>::operator=( const ComplexBuildVectorExpression<RealType<ValueType> >& exp )
{
    buildComplex( exp.getRealArg(), exp.getImagArg() );
    return *this;
}

template<typename ValueType>
Vector<ValueType>& Vector<ValueType>::operator+=( const Vector<ValueType>& other )
{
    setVector( other, common::BinaryOp::ADD );
    return *this;
}

template<typename ValueType>
template<typename OtherValueType>
Vector<ValueType>& Vector<ValueType>::operator+=( const CastVectorExpression<ValueType, OtherValueType>& other )
{
    setVector( other.getArg(), common::BinaryOp::ADD );
    return *this;
}

template<typename ValueType>
Vector<ValueType>& Vector<ValueType>::operator-=( const Vector<ValueType>& other )
{
    setVector( other, common::BinaryOp::SUB );
    return *this;
}

template<typename ValueType>
template<typename OtherValueType>
Vector<ValueType>& Vector<ValueType>::operator-=( const CastVectorExpression<ValueType, OtherValueType>& other )
{
    setVector( other.getArg(), common::BinaryOp::SUB );
    return *this;
}

template<typename ValueType>
Vector<ValueType>& Vector<ValueType>::operator*=( const Vector<ValueType>& other )
{
    setVector( other, common::BinaryOp::MULT );
    return *this;
}

template<typename ValueType>
template<typename OtherValueType>
Vector<ValueType>& Vector<ValueType>::operator*=( const CastVectorExpression<ValueType, OtherValueType>& other )
{
    setVector( other.getArg(), common::BinaryOp::MULT );
    return *this;
}

template<typename ValueType>
Vector<ValueType>& Vector<ValueType>::operator/=( const Vector<ValueType>& other )
{
    setVector( other, common::BinaryOp::DIVIDE );
    return *this;
}

template<typename ValueType>
template<typename OtherValueType>
Vector<ValueType>& Vector<ValueType>::operator/=( const CastVectorExpression<ValueType, OtherValueType>& other )
{
    setVector( other.getArg(), common::BinaryOp::DIVIDE );
    return *this;
}

/* ------------------------------------------------------------------------- */
/*  Inline translation:  operator=  scalar/vector binop scalar/vector        */
/* ------------------------------------------------------------------------- */

template<typename ValueType>
template<common::BinaryOp op>
Vector<ValueType>& Vector<ValueType>::operator=( const Expression<intern::Scalar, Vector<ValueType>, op>& exp )
{
    intern::Scalar alpha = exp.getArg1();
    this->binaryOpScalar( exp.getArg2(), alpha.getValue<ValueType>(), op, true );  // swap 
    return *this;
}

template<typename ValueType>
template<common::BinaryOp op>
Vector<ValueType>& Vector<ValueType>::operator=( const Expression<Vector<ValueType>, intern::Scalar, op>& exp )
{
    intern::Scalar alpha = exp.getArg2();
    this->binaryOpScalar( exp.getArg1(), alpha.getValue<ValueType>(), op, false );  // no swap
    return *this;
}

template<typename ValueType>
template<common::BinaryOp op>
Vector<ValueType>& Vector<ValueType>::operator=( const Expression<Vector<ValueType>, Vector<ValueType>, op>& exp )
{
    this->binaryOp( exp.getArg1(), op, exp.getArg2() );  
    return *this;
}

/* ------------------------------------------------------------------------- */
/*  Inline translation:  operator=  unaryFunction( x )                       */
/* ------------------------------------------------------------------------- */

template<typename ValueType>
Vector<ValueType>& Vector<ValueType>::operator=( const UnaryVectorExpression<ValueType>& unaryVectorExp )
{
    unaryOp( unaryVectorExp.getArg(), unaryVectorExp.getOp() );
    return* this;
}

template<typename ValueType>
void Vector<ValueType>::unaryOpInPlace( const common::UnaryOp op )
{
    unaryOp( *this, op );
}

template<typename ValueType>
void Vector<ValueType>::binaryOp( const ValueType& alpha, const common::BinaryOp op, const Vector<ValueType>& x )
{
    binaryOpScalar( x, alpha, op, true );    // swap order of the arguments
}

template<typename ValueType>
void Vector<ValueType>::binaryOp( const Vector<ValueType>& x, const common::BinaryOp op, const ValueType& alpha ) 
{
    binaryOpScalar( x, alpha, op, false );   // keep order of the arguments
}

} /* end namespace lama */

} /* end namespace scai */
