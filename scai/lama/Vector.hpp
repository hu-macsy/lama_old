/**
 * @file Vector.hpp
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
 * @brief Definition of an abstract class for distributed vectors.
 * @author Jiri Kraus
 * @date 22.02.2011
 */
#pragma once

// for dll_import
#include <scai/common/config.hpp>

// base classes
#include <scai/dmemo/Distributed.hpp>

// local library
#include <scai/lama/expression/Expression.hpp>

#include <scai/lama/Scalar.hpp>
#include <scai/lama/io/FileIO.hpp>

// others
#include <scai/hmemo.hpp>

#include <scai/logging.hpp>

#include <scai/common/Factory.hpp>
#include <scai/common/ScalarType.hpp>
#include <scai/common/SCAITypes.hpp>

#include <utility>

namespace scai
{

namespace lama
{

class Matrix;

/** Pointer class for a vector, always use of a shared pointer. */

typedef common::shared_ptr<class Vector> VectorPtr;

/** Help class as forward declaration of enum types belonging to class Vector. */

struct _Vector
{

    /**
     * @brief VectorFormat describes if a vector is dense or sparse.
     */
    typedef enum
    {
        DENSE,      //!< vector format for a dense vector
        SPARSE,     //!< vector format for a sparse vector, not supported yet
        UNDEFINED   //!< for convenience, always the last entry, stands also for number of entries
    } VectorFormat;

    static COMMON_DLL_IMPORTEXPORT const char* kind2Str( const VectorFormat vectorKind );

    static COMMON_DLL_IMPORTEXPORT VectorFormat str2Kind( const char* str );

};  // struct _Vector

/** @brief Output operator<< for VectorFormat prints meaningful names instead of int values */

COMMON_DLL_IMPORTEXPORT std::ostream& operator<<( std::ostream& stream, const _Vector::VectorFormat& kind );

/** Type definition for the key type used for the Vector factory.
 *
 *  The key for vector create is a pair of vector format and the value type.
 */

typedef std::pair<_Vector::VectorFormat, common::scalar::ScalarType> VectorCreateKeyType;

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
    public dmemo::Distributed,
    public _Vector

{
public:

    /** @brief More convenient use of the create routine of factory that avoids use of CreateKeyType.
     */
    static Vector* getVector( const VectorFormat format, const common::scalar::ScalarType valueType );

    /** @brief Create a dense vector of a certain value type and a given distribution.
     *
     *  This method keeps compatibility with an older method that did know which vectors were supported.
     */
    static Vector* getDenseVector( const common::scalar::ScalarType valueType, dmemo::DistributionPtr distribution );

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

    /** this = x * y */

    Vector& operator=( const Expression_VV );

    /** this = alpha * x * y */

    Vector& operator=( const Expression_SVV );

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
     * @brief Multiplies the passed value with all elements of this.
     *
     * @param[in] other   the vector to multiply to do the multiplication per element
     * @return            a reference to this.
     */
    Vector& operator*=( const Vector& other );

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
    virtual void buildValues( hmemo::_HArray& values ) const = 0;

    /**
     * @brief Sets the local values of a vector by an array.
     *
     * @param[out] values    is the array with local vector values.
     *
     * Note: A conversion operator must be available for values.getValueType() to
     *       the type of this vector.
     */
    virtual void setValues( const hmemo::_HArray& values ) = 0;

    /**
     * This method initializes a (distributed) vector with random numbers. 
     * 
     * @param[in] distribution specifies the distribution of the vector
     * @param[in] fillRate for the number of non-zeros
     */
    virtual void setRandom( dmemo::DistributionPtr distribution, const float fillRate = 1.0 ) = 0;

    /**
     * This method sets a vector by reading its values from one or multiple files.
     *
     * @param[in] fileName      the filename to read from
     * @param[in] distribution  optional, if set it is the distribution of the vector 
     *
     *   \code
     *      DenseVector<double> vector;
     *      vector.readFromFile( "vector.mtx" )                    ! vector only on processor 0
     *      vector.readFromFile( "vector_%r.mtx" )                 ! general block distributed vector, each processor reads it own file
     *      vector.readFromFile( "vector.mtx", rowDist )           ! each processor gets its local part of the vector in one file
     *      vector.readFromFile( "vector_%r.mtx", rowDist )        ! read a partitioned vector with the given distribution
     *   \endcode
     */
    void readFromFile( const std::string& fileName, dmemo::DistributionPtr distribution = dmemo::DistributionPtr() );

    /**
     *  This method sets a vector a reading its values from one or multiple files and also the distribution from a file
     *
     * @param[in] vectorFileName the single or partitioned filename to read from
     * @param[in] distributionFileName the single or partitioned filename with the row distribution of the vector
     *
     *   \code
     *      CSRSparseMatrix<double> vector;
     *      vector.readFromFile( "vector.mtx", "owners.mtx" )
     *      vector.readFromFile( "vector_%r.mtx", "owners.mtx" )
     *      vector.readFromFile( "vector.mtx", "rows%r.mtx" )
     *      vector.readFromFile( "vector_%r.mtx", "rows%r.mtx" )
     *   \endcode
     */
    void readFromFile( const std::string& vectorFileName, const std::string& distributionFileName );

    /**
     * @brief write the vector to an output file
     *
     * @param[in] fileName is the name of the output file (suffix must be added according to the file type)
     * @param[in] fileType format of the output file ("frv" for SAMG, "mtx" for MatrixMarket), default is to decide by suffix
     * @param[in] dataType representation type for output values, if set it overrides IO settings
     * @param[in] fileMode can be BINARY or FORMATTED, DEFAULT_MODE keeps default/environment settings
     */
    void writeToFile(
        const std::string& fileName,
        const std::string& fileType = "",
        const common::scalar::ScalarType dataType = common::scalar::UNKNOWN,
        const FileIO::FileMode fileMode = FileIO::DEFAULT_MODE  ) const;

    /**
     * @brief get a vector with all local values
     */
    virtual const hmemo::_HArray& getLocalValues() const = 0;

    /**
     * @brief Queries the value type of the vector elements, e.g. DOUBLE or FLOAT.
     */
    virtual common::scalar::ScalarType getValueType() const = 0;

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
     *  Method to create a new vector of the same kind and same type
     *   
     *  /code
     *    const Vector& old = ...
     *    Vector* new = Vector::create( old.getCreateValue() );
     *  /endcode
     *
     *  This routine is very important to write code that can deal with arbitrary types
     *  but does not have a template param for the value type.
     */

    virtual VectorCreateKeyType getCreateValue() const = 0;

    /**
     *  @brief Creates a new Vector of the same kind and value type, and same context
     *
     *  /code
     *    const Vector& old = ...
     *    ....
     *    Vector* new = old.newVector();
     *
     *    // is same as
     *
     *    Vector* new = Vector::create( old.getCreateValue() );
     *    new->setContextPtr( new.getContextPtr() ); 
     *  /endcode
     *
     *  The new vector is a zero vector, not allocated, not initialized.
     */
    virtual Vector* newVector() const = 0;

    /**
     *  @brief copy is a virtual call of the copy constructor of the derived classes
     *
     *  /code
     *    const Vector& old = ...
     *    Vector* new = old.cooy()
     *
     *    // is same as
     *
     *    Vector* new = Vector::create( old.getCreateValue() );
     *    *new = old;
     *  /endcode
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
     * Swap is only possible if both vectors are of the same format (DENSE) and
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
    virtual void assign( const hmemo::_HArray& localValues, dmemo::DistributionPtr distribution ) = 0;

    /**
     *  Builds an array with local values of a distributed vector.
     *
     *  @param[out] localValues   will be an array that contains local values of the vector
     *
     *  For different value types, implicit format conversion will be done.
     *  A sparse vector should generate an array with all values.
     */
    virtual void buildLocalValues( hmemo::_HArray& localValues ) const = 0;

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

    virtual void assign( const Expression_SVV& expression ) = 0;

    /**
     * @brief Returns the dot product of this and other.
     *
     * @param[in] other   the vector to calculate the dot product with.
     * @return            the dot product of this and other
     */
    virtual Scalar dotProduct( const Vector& other ) const = 0;

    /**
     *  @brief Scale a Vector with another Vector.
     *
     *  @param[in] other   the other vector to scale this with
     *  @return            reference to the scaled vector
     */
    virtual Vector& scale( const Vector& other ) = 0;

    /**
     * @brief Starts a prefetch to make this valid at the passed context.
     *
     * @param[in] context specifies the location to make this vector valid at
     */
    virtual void prefetch( const hmemo::ContextPtr context ) const = 0;

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
    void setContextPtr( hmemo::ContextPtr context );

    /**
     * @brief Getter function for the context (pointer) of a vector.
     */
    inline hmemo::ContextPtr getContextPtr() const;

    /**
     * @brief Returns the global memory that is allocated to hold this vector.
     * For a distributed vector all partitions are summed together.
     *
     * @return the memory consumption of this vector.
     */
    virtual size_t getMemoryUsage() const = 0;

    /**
     *  @brief Allocates or reallocates this vector for a given distribution.
     *
     *  All elements of the vector are undefined after this operation.
     */
    virtual void allocate( dmemo::DistributionPtr distributionPtr ) = 0;

    /**
     * @brief Redistributes this vector to the new passed distribution.
     *
     * @param[in] distribution   the new distribution for this vector.
     *
     * The global vector itself remains unchanged; only local parts
     * can be different now.
     */
    virtual void redistribute( dmemo::DistributionPtr distribution ) = 0;

    /**
     *  Build the conjugate vector in place.
     */
    virtual void conj() = 0;

    /**
     *  Calculates the exponentional function of the vector elements in place.
     */
    virtual void exp() = 0;

protected:

    /**
     *  Constructor of Vector for derived classes by size and/or context
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

    hmemo::ContextPtr mContext; //!< decides about location of vector operations

    SCAI_LOG_DECL_STATIC_LOGGER( logger )

private:

    /** write only the local data to a file, no communication here */

    void writeLocalToFile(
        const std::string& fileName,
        const std::string& fileType,
        const common::scalar::ScalarType dataType,
        const FileIO::FileMode fileMode ) const;

    /** write the whole vector into a single file, can imply redistribution */

    void writeToSingleFile(
        const std::string& fileName,
        const std::string& fileType,
        const common::scalar::ScalarType dataType,
        const FileIO::FileMode fileMode ) const;

    /** same as writeLocalToFile but also communication for error handling */

    void writeToPartitionedFile(
        const std::string& fileName,
        const std::string& fileType,
        const common::scalar::ScalarType dataType,
        const FileIO::FileMode fileMode ) const;

    void readFromSingleFile( const std::string& fileName );

    void readFromPartitionedFile( const std::string& myPartitionFileName, dmemo::DistributionPtr dist );
};

IndexType Vector::size() const
{
    return getDistributionPtr()->getGlobalSize();
}

hmemo::ContextPtr Vector::getContextPtr() const
{
    return mContext;
}

} /* end namespace lama */

} /* end namespace scai */
