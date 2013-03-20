/**
 * @file Vector.hpp
 *
 * @license
 * Copyright (c) 2011
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
 * @brief Vector.hpp
 * @author Jiri Kraus
 * @date 22.02.2011
 * $Id$
 */
#ifndef LAMA_VECTOR_HPP_
#define LAMA_VECTOR_HPP_

// for dll_import
#include <lama/config.hpp>

// base classes
#include <lama/Distributed.hpp>

// others
#include <lama/expression/Expression.hpp>
#include <lama/ContextFactory.hpp>

#include <lama/LAMATypes.hpp>
#include <lama/Scalar.hpp>
#include <lama/Context.hpp>

// logging
#include <logging/logging.hpp>

namespace lama
{

class Matrix;

/** Pointer class for a vector, always use of a shared pointer. */

typedef boost::shared_ptr<class Vector> VectorPtr;

/**
 * @brief The class Vector is a abstract type that represents a distributed 1D real or complex vector.
 *
 * Vector is one of the LAMA base types and should be used in all situations where it is not necessary to access a
 * single element or to create a new Vector.
 *
 * As this class is an abstract class, all constructors are protected.
 *
 * The following methods must be implemented by derived classes:
 *
 *  - buildValues to get all local elements on one processor
 *  - getLocalValue( i ) the the i-th local element
 *
 * This base class can be used to define dense and sparse vectors of
 * any type.
 */
class LAMA_DLL_IMPORTEXPORT Vector: public Distributed
{
public:

    static std::auto_ptr<Vector> createVector( const Scalar::ScalarType valueType, DistributionPtr distribution );

    /**
     * @brief ExpressionMemberType is the type that is used the template Expression to store a Vector.
     */
    typedef const Vector& ExpressionMemberType;

    /**
     * @brief Releases all allocated resources.
     */
    virtual ~Vector();

    /**
     * @brief The assignment operator assigns the result of the passed expression
     *        to this.
     *
     * The assignment operator assigns the result of the passed expression to
     * this, if necessary new memory will be allocated. The Vector will hold
     * the result of the Matrix Vector multiplication represented by
     * expression.
     *
     * @param[in] expression the input expression.
     * @return               a reference to this.
     * @throws               Exceptions thrown by the Allocator
     */
    Vector& operator=( const Expression<Matrix,Vector,Times>& expression );

    Vector& operator=( const Expression<Scalar,Expression<Matrix,Vector,Times>,Times>& expression );

    Vector& operator=(
        const Expression<Expression<Scalar,Vector,Times>,Expression<Scalar,Vector,Times>,Plus>& expression );

    Vector& operator=(
        const Expression<Expression<Scalar,Expression<Matrix,Vector,Times>,Times>,Expression<Scalar,Vector,Times>,Plus>& expression );

    Vector& operator=( const Expression<Scalar,Vector,Times>& expression );

    Vector& operator=( const Expression<Vector,Vector,Plus>& expression );

    Vector& operator=( const Vector& other );

    Vector& operator*=( const Scalar value );

    Vector& operator+=( const Vector& other );

    Vector& operator+=( const Expression<Scalar,Vector,Times>& expression );

    Vector& operator-=( const Vector& other );

    /**
     * @brief Assigns the passed value to all elements of this.
     *
     * @param[in] value the value to assign to all elements of this.
     * @return          a reference to this.
     */
    Vector& operator=( const Scalar value );

    /**
     * @brief returns a copy of the value at the passed global index.
     *
     * @param[in] i the global index to get the value at.
     * @return       a copy of the value at the passed global position.
     *
     * As this operator requires communication ins SPMD mode it can be very inefficient in some situations.
     */
    const Scalar operator()( const IndexType i ) const;

    /**
     * @brief returns the dot product of this and other
     *
     * @param[in] other the vector do calculate the dot product with.
     * @return          the dot product of this and other
     */
    Scalar operator*( const Vector& other ) const;

    /**
     * @brief Build an array with local values of the vector.
     *
     * @param[in,out] values LAMA array that will be filled with the local values.
     *
     * Only the type of the LAMA array is used as input arg to determine the value type.
     */
    virtual void buildValues( _LAMAArray& values ) const = 0;

    /**
     * @brief Set the local values of a vector by an array.
     *
     * @param[out] values is the array with local vector values.
     *
     * Note: A conversion operator must be available for values.getValueType() to
     *       the type of this vector.
     */
    virtual void setValues( const _LAMAArray& values ) = 0;

    /**
     * @brief Query the value type of the vector elements, e.g. DOUBLE or FLOAT.
     */
    virtual Scalar::ScalarType getValueType() const = 0;

    /**
     * @brief returns a copy of the value at the passed global index.
     *
     * @param[in] globalIndex   the global index to get the value at.
     * @return                  a copy of the value at the passed global position.
     *
     * As this operation requires communication in SPMD mode it can be very inefficient in some situations.
     */
    virtual Scalar getValue( IndexType globalIndex ) const =0;

    /**
     * @brief returns the global minimum value of this.
     *
     * @return the global minimum value of this vector.
     */
    virtual Scalar min() const =0;

    /**
     * @brief returns the global maximum value of this.
     *
     * @return the global maximum value of this vector.
     */
    virtual Scalar max() const =0;

    /**
     * @brief returns the L1 norm of this.
     *
     * @return the L1 norm of this.
     *
     * l1Norm computes the sum of the absolute values of this.
     */
    virtual Scalar l1Norm() const =0;

    /**
     * @brief returns the L2 norm of this.
     *
     * @return the L2 norm of this.
     *
     * l2Norm computes the sum of the absolute values of this.
     */
    virtual Scalar l2Norm() const =0;

    /**
     * @brief returns the max norm of this.
     *
     * @return the max norm of this.
     *
     * maxNorm computes the value of this with the largest magnitude.
     */
    virtual Scalar maxNorm() const =0;

    //TODO: This is an uninitialized create, do we need a initialized create?
    //      What is a good interface for an initialized create?
    /**
     * @brief Create is a virtual constructor, which creates a new Vector with the same concrete class, size and distribuiton than this.
     *
     * @return                  a pointer to the new Vector, caller has the owner ship.
     */
    virtual std::auto_ptr<Vector> create() const =0;

    /**
     * @brief Create is a virtual constructor, which creates a new Vector with the same concrete class than this.
     *
     * @param[in] distribution  the distribution to use for the new Vector.
     * @return                  a pointer to the new Vector, caller has the owner ship.
     */
    virtual std::auto_ptr<Vector> create( DistributionPtr distribution ) const =0;

    /**
     * @brief Returns the size of the vector.
     *
     * @return  the size of this vector.
     */
    inline IndexType size() const;

    /**
     * @brief Swap the content of of this vector with another vector.
     *
     * @param[in,out] other the Vector to swap the contents with.
     *
     * Swap is only possible if both vectors are of the same kind (DENSE) and
     * have the same value type.
     */
    virtual void swap( Vector& other ) = 0;

    /**
     * @brief Write some information about this to the passed stream.
     *
     * @param[out] stream   the stream to write to.
     */
    virtual void writeAt( std::ostream& stream ) const;

    /**
     *  @brief Assign an arbitrary vector to this vector.
     */

    virtual void assign( const Vector& other ) =0;

    /**
     *  Assignment to vector by local values and distribution.
     */

    virtual void assign( const _LAMAArray& localValues, DistributionPtr distribution ) = 0;

    /**
     *  Build an array with local values of a distributed vector.
     *
     *  @param[out] localValues will be an array that contains local values of the vector
     *
     *  For different value types, implicit format conversion will be done.
     *  A sparse vector should generate an array with all values.
     */

    virtual void buildLocalValues( _LAMAArray& localValues ) const = 0;

    /**
     * @brief Assigns the passed value to all elements of this.
     *
     * @param[in] value the value to assign to all elements of this.
     */
    virtual void assign( const Scalar value ) =0;

    /**
     * @brief Assignment of a 'full' vector expression.
     */
    virtual void assign(
        const Expression<Expression<Scalar,Vector,Times>,Expression<Scalar,Vector,Times>,Plus>& expression ) =0;

    virtual Scalar dotProduct( const Vector& other ) const =0;

    /**
     * @brief Starts a prefetch to make this valid at the passed context.
     *
     * @param[in] context specifies the location to make this vector valid at
     */
    virtual void prefetch( const ContextPtr context ) const =0;

    /**
     * @brief wait for a possibly running prefetch.
     */
    virtual void wait() const =0;

    /**
     * @brief This method inverts all elements of the vector and is completely local.
     */

    virtual void invert() =0;

    /**
     * @brief Set the 'preferred' context where data resides and computations are done
     */
    void setContext( ContextPtr location );

    /**
     * @brief Getter for the context (pointer) of a vector.
     */
    inline ContextPtr getContext() const;

    /**
     * @brief getMemoryUsage returns the global memory that is allocated to hold this vector.
     *
     * getMemoryUsage returns the global memory that is allocated to hold this vector. For a distributed vector
     * all partitions are summed together.
     *
     * @return the memory consumption of this vector.
     */
    virtual size_t getMemoryUsage() const =0;

    /**
     *  @brief Allocate this vector for a given distribution.
     *
     *  All elements of the vector are undefined after this operation.
     *  Elements can be set e.g. with
     */
    void resize( DistributionPtr distributionPtr );

    /**
     * @brief redistributes this vector to the new passed distribution.
     *
     * @param[in] distribution  the new distribution for this vector.
     *
     * The global vector itself remains unchanged; only local parts
     * can be now different.
     */
    virtual void redistribute( DistributionPtr distribution ) =0;

protected:

    /**
     *  Constructor of Vector for derived classes by size and/or context
     */

    explicit Vector( const IndexType size = 0, ContextPtr context = ContextFactory::getContext( Context::Host ) );

    /**
     * @brief Constructor of Vector for derived classes by distribution
     *
     * @param[in] distribution  the distribution to use for the new Vector.
     * @param[in] context is optional, will be Host context
     */
    explicit Vector( DistributionPtr distribution, ContextPtr context = ContextFactory::getContext( Context::Host ) );

    /**
     * @brief Creates a copy of the passed Vector.
     *
     * @param[in] other the Vector to take a copy from.
     *
     * Inherits size/distribution and the context of the passed vector.
     */
    Vector( const Vector& other );

    /**
     *  @brief Swap member variables of Vector class.
     */

    void swapVector( Vector& other );

    virtual void resizeImpl() =0;

    ContextPtr mContext; //!< decides about location of vector operations

    LAMA_LOG_DECL_STATIC_LOGGER( logger );
};

IndexType Vector::size() const
{
    return getDistributionPtr()->getGlobalSize();
}

ContextPtr Vector::getContext() const
{
    return mContext;
}

}

#endif // LAMA_VECTOR_HPP_
