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
#pragma once

// for dll_import
#include <common/config.hpp>

// base classes
#include <lama/Distributed.hpp>

// others
#include <lama/expression/Expression.hpp>

#include <memory.hpp>

#include <lama/LAMATypes.hpp>
#include <lama/Scalar.hpp>
#include <lama/io/FileType.hpp>

// logging
#include <logging.hpp>

#include <common/Factory.hpp>
#include <common/ScalarType.hpp>

// Not recommened, but here for convenience

using memory::ContextPtr;

namespace lama
{

class Matrix;

/** Pointer class for a vector, always use of a shared pointer. */

typedef common::shared_ptr<class Vector> VectorPtr;

class Vector;

/**
 * @brief VectorKind describes if a vector is dense or sparse.
 */
typedef enum
{
    DENSE, //!< vector kind for a dense vector
    SPARSE //!< vector kind for a sparse vector, not supported yet

} VectorKind;

/** For convenience: add the key type used for the Vector factory. */

typedef std::pair<VectorKind, common::ScalarType> VectorCreateKeyType;

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
class COMMON_DLL_IMPORTEXPORT Vector: 

     public common::Factory<VectorCreateKeyType, Vector*>,
     public Distributed

{
public:

    /**
     * @brief Vector factory to get a vector of a certain kind and a certain type
     *
     * @param[in] kind is either DENSE or SPARSE
     * @param[in] valueType specifies the value type as the elements, e.g. FLOAT, DOUBLE
     *
     * This factory operation allows to create a vector at runtime of any format or any type.
     * Internally, all vector classes must register their create operation.
     */
    static Vector* getVector( const VectorKind kind, const common::ScalarType valueType );

    /** @brief Create a dense vector of a certain value type and a given distribution.
     *
     *  This method keeps compatibility with an older method that did know which vectors were supported.
     */
    static Vector* createVector( const common::ScalarType valueType, DistributionPtr distribution );

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
    Vector& operator=( const Expression_MV& expression );

    Vector& operator=( const Expression_VM& expression );

    /** this = alpha * A * x */

    Vector& operator=( const Expression_SMV& expression );

    /** this = alpha * x * A */

    Vector& operator=( const Expression_SVM& expression );

    /** this = alpha * x + beta * y */

    Vector& operator=( const Expression_SV_SV& expression );

    /** this = alpha * A * x + beta * y */

    Vector& operator=( const Expression_SMV_SV& expression );

    /** this = alpha * x * A + beta * y */

    Vector& operator=( const Expression_SVM_SV& expression );

    /** this = alpha * x */

    Vector& operator=( const Expression_SV& expression );

    /** this +=  alpha * A * x */

    Vector& operator+=( const Expression_SMV& expression );

    Vector& operator+=( const Expression_SVM& expression );

    /** this +=  alpha * x */

    Vector& operator+=( const Expression_SV& expression );

    /** this -=  alpha * A * x */

    Vector& operator-=( const Expression_SMV& expression );

    Vector& operator-=( const Expression_SVM& expression );

    /** this -=  alpha * x */

    Vector& operator-=( const Expression_SV& expression );

    /**
     * @brief Assigns the values of other to the elements of this.
     *
     * @param[in] other   the vector to get values from.
     * @return            a reference to this.
     */
    Vector& operator=( const Vector& other );

    /**
     * @brief Multiplies the passed value with all elements of this.
     *
     * @param[in] value   the value to multiply all elements of this with.
     * @return            a reference to this.
     */
    Vector& operator*=( const Scalar value );

    /**
     * @brief Divides the passed value with all elements of this.
     *
     * @param[in] value   the value to divide all elements of this with.
     * @return            a reference to this.
     */
    Vector& operator/=( const Scalar value );

    /**
     * @brief Returns the addition of this and other.
     *
     * @param[in] other the vector to do the addition with.
     * @return          a reference to this.
     */
    Vector& operator+=( const Vector& other );

    /**
     * @brief Returns the subtraction of this and other.
     *
     * @param[in] other the vector to do the subtraction with.
     * @return          a reference to this.
     */
    Vector& operator-=( const Vector& other );

    /**
     * @brief Assigns the passed value to all elements of this.
     *
     * @param[in] value   the value to assign to all elements of this.
     * @return            a reference to this.
     */
    Vector& operator=( const Scalar value );

    /**
     * @brief Returns a copy of the value at the passed global index.
     *
     * @param[in] i    the global index to get the value at.
     * @return         a copy of the value at the passed global position.
     *
     * As this operator requires communication ins SPMD mode it can be very inefficient in some situations.
     */
    const Scalar operator()( const IndexType i ) const;

    /**
     * @brief Builds an array with local values of the vector.
     *
     * @param[in,out] values   LAMA array that will be filled with the local values.
     *
     * Only the type of the LAMA array is used as input arg to determine the value type.
     */
    virtual void buildValues( memory::ContextArray& values ) const = 0;

    /**
     * @brief Sets the local values of a vector by an array.
     *
     * @param[out] values    is the array with local vector values.
     *
     * Note: A conversion operator must be available for values.getValueType() to
     *       the type of this vector.
     */
    virtual void setValues( const memory::ContextArray& values ) = 0;

    /**
     * @brief Assign this vector with values stored the file with the given filename.
     *
     * @param[in] filename  the name of the file to be read containing vector data.
     *
     * The implementation of this method in derived classes can make its own
     * decision what the distribtion of this vector will be.
     */
    virtual void readFromFile( const std::string& filename ) = 0;

    virtual void writeToFile(
        const std::string& fileBaseName,
        const File::FileType fileType = File::XDR,
        const File::DataType dataType = File::DOUBLE ) const = 0;

    /**
     * @brief get a vector with all local values
     */
    virtual const memory::ContextArray& getLocalValues() const = 0;

    /**
     * @brief Queries the value type of the vector elements, e.g. DOUBLE or FLOAT.
     */
    virtual common::ScalarType getValueType() const = 0;

    /**
     * @brief Returns a copy of the value at the passed global index.
     *
     * @param[in] globalIndex   the global index to get the value at.
     * @return                  a copy of the value at the passed global position.
     *
     * As this operation requires communication in SPMD mode it can be very inefficient in some situations.
     */
    virtual Scalar getValue( IndexType globalIndex ) const = 0;

    /**
     * @brief Returns the global minimum value of this.
     *
     * @return   the global minimum value of this vector.
     */
    virtual Scalar min() const = 0;

    /**
     * @brief Returns the global maximum value of this.
     *
     * @return the global maximum value of this vector.
     */
    virtual Scalar max() const = 0;

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
     * @brief Returns the max norm of this.
     *
     * @return the max norm of this.
     *
     * maxNorm computes the value of this with the largest magnitude.
     */
    virtual Scalar maxNorm() const = 0;

    /**
     * @brief clone is a virtual call of the default constructor of the derived classes
     *
     * @return a pointer to the new Vector, caller takes the ownership.
     */
    virtual Vector* clone() const = 0;

    /**
     * @brief Create is a virtual constructor, which clones a new Vector with the same concrete class as this.
     *
     * @param[in] distribution  the distribution to use for the new Vector.
     * @return                  a pointer to the new Vector, caller has the owner ship.
     */
    virtual Vector* clone( DistributionPtr distribution ) const = 0;

    /**
     *  @brief copy is a virtual call of the copy constructor of the derived classes
     */
    virtual Vector* copy() const = 0;

    /**
     * @brief Returns the size of the vector.
     *
     * @return  the size of this vector.
     */
    inline IndexType size() const;

    /**
     * @brief Swaps the content of this vector with another vector.
     *
     * @param[in,out] other   the Vector to swap the contents with.
     *
     * Swap is only possible if both vectors are of the same kind (DENSE) and
     * have the same value type.
     */
    virtual void swap( Vector& other ) = 0;

    virtual void writeAt( std::ostream& stream ) const;

    /**
     *  @brief Assigns an arbitrary vector to this vector.
     */
    virtual void assign( const Vector& other ) = 0;

    /**
     *  Assignment to vector by local values and distribution.
     */
    virtual void assign( const memory::ContextArray& localValues, DistributionPtr distribution ) = 0;

    /**
     *  Builds an array with local values of a distributed vector.
     *
     *  @param[out] localValues   will be an array that contains local values of the vector
     *
     *  For different value types, implicit format conversion will be done.
     *  A sparse vector should generate an array with all values.
     */
    virtual void buildLocalValues( memory::ContextArray& localValues ) const = 0;

    /**
     * @brief Assigns the passed value to all elements of this.
     *
     * @param[in] value   the value to assign to all elements of this.
     */
    virtual void assign( const Scalar value ) = 0;

    /**
     * @brief Assignment of a 'full' vector expression.
     */
    virtual void assign( const Expression_SV_SV& expression ) = 0;

    /**
     * @brief Returns the dot product of this and other.
     *
     * @param[in] other   the vector to calculate the dot product with.
     * @return            the dot product of this and other
     */
    virtual Scalar dotProduct( const Vector& other ) const = 0;

    /**
     * @brief Starts a prefetch to make this valid at the passed context.
     *
     * @param[in] context specifies the location to make this vector valid at
     */
    virtual void prefetch( const ContextPtr context ) const = 0;

    /**
     * @brief Starts a prefetch to make data valid at the context of the vector.
     *
     */

    void prefetch() const;

    /**
     * @brief Waits for a possibly running prefetch.
     */
    virtual void wait() const = 0;

    /**
     * @brief This method inverts all elements of the vector and is completely local.
     */
    virtual void invert() = 0;

    /**
     * @brief Sets the 'preferred' context where data resides and computations are done.
     */
    void setContext( ContextPtr location );

    /**
     * @brief Getter function for the context (pointer) of a vector.
     */
    inline ContextPtr getContext() const;

    /**
     * @brief Returns the global memory that is allocated to hold this vector.
     * For a distributed vector all partitions are summed together.
     *
     * @return the memory consumption of this vector.
     */
    virtual size_t getMemoryUsage() const = 0;

    /**
     *  @brief Allocates this vector for a given distribution.
     *
     *  All elements of the vector are undefined after this operation.
     *  Elements can be set e.g. with
     */
    void resize( DistributionPtr distributionPtr );

    /**
     * @brief Redistributes this vector to the new passed distribution.
     *
     * @param[in] distribution   the new distribution for this vector.
     *
     * The global vector itself remains unchanged; only local parts
     * can be different now.
     */
    virtual void redistribute( DistributionPtr distribution ) = 0;

protected:

    /**
     *  Constructor of Vector for derived classes by size and/or context
     */
    explicit Vector( const IndexType size = 0, ContextPtr context = ContextPtr() );

    /**
     * @brief Constructor of Vector for derived classes by distribution
     *
     * @param[in] distribution  the distribution to use for the new Vector.
     * @param[in] context       is optional, will be Host context.
     */
    explicit Vector( DistributionPtr distribution, ContextPtr context = ContextPtr() );

    /**
     * @brief Creates a copy of the passed Vector.
     *
     * @param[in] other   the Vector to take a copy from.
     *
     * Inherits size/distribution, context and content of the passed vector.
     */
    Vector( const Vector& other );

    /**
     *  @brief Swaps member variables of Vector class.
     */
    void swapVector( Vector& other );

    /**
     *  @brief TODO[doxy] Complete Description.
     */
    virtual void resizeImpl() = 0;

    ContextPtr mContext; //!< decides about location of vector operations

    LAMA_LOG_DECL_STATIC_LOGGER( logger )
};

IndexType Vector::size() const
{
    return getDistributionPtr()->getGlobalSize();
}

ContextPtr Vector::getContext() const
{
    return mContext;
}

}  // namespace

/** @brief  stream output for key values of creator  */

inline std::ostream& operator<<( std::ostream& stream, const lama::VectorCreateKeyType& key )
{
    stream << "<" << key.first << ", " << key.second << ">";
    return stream;
}

