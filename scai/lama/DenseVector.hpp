/**
 * @file DenseVector.hpp
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
 * @brief Definition of template class that stands for a dense vector of a certain type.
 * @author Jiri Kraus
 * @date 22.02.2011
 */

#pragma once

// for dll_import
#include <scai/common/config.hpp>

// base classes
#include <scai/lama/Vector.hpp>

// internal scai libraries
#include <scai/utilskernel/LArray.hpp>
#include <scai/dmemo/Distribution.hpp>
#include <scai/dmemo/Halo.hpp>
#include <scai/hmemo.hpp>

#include <scai/tasking/SyncToken.hpp>

#include <scai/common/macros/throw.hpp>

// std
#include <fstream>

namespace scai
{

namespace lama
{

/**
 * @brief The template DenseVector represents a distributed 1D Vector with elements of type ValueType.
 *
 * @tparam ValueType the value type for the vector values.
 */
template<typename ValueType>
class COMMON_DLL_IMPORTEXPORT DenseVector:

    public Vector,

    public Vector::Register<DenseVector<ValueType> >    // register at factory

{
public:

    /** Default constructor, creates empty (not initilized) vector, that is replicated (without distribution) */

    DenseVector();

    /**
     * @brief creates a not initialized distributed DenseVector of the passed global size.
     *
     * @param[in] context  the context to use for the new vector.
     */
    explicit DenseVector( hmemo::ContextPtr context );

    /**
     * @brief creates a not initialized distributed DenseVector of the passed global size.
     *
     * @param[in] distribution  the distribution to use for the new vector.
     */
    explicit DenseVector( dmemo::DistributionPtr distribution );

    /**
     * @brief creates a not initialized distributed DenseVector of the passed global size.
     *
     * @param[in] distribution  the distribution to use for the new vector.
     * @param[in] context  the context to use for the new vector.
     */
    explicit DenseVector ( dmemo::DistributionPtr distribution, hmemo::ContextPtr context );

    /**
     * @brief creates a replicated DenseVector of the passed size initialized to the passed value.
     *
     * @param[in] size  the size of the new DenseVector.
     * @param[in] value the value to assign to all elements of the new DenseVector.
     * @param[in] context   specifies optionally the context where dense vector should reside
     */
    DenseVector( const IndexType size, const ValueType value, hmemo::ContextPtr context = hmemo::ContextPtr() );

    /**
     * @brief creates a distributed DenseVector of the passed global size initialized to the passed value.
     *
     * @param[in] distribution  the distribution to use for the new vector.
     * @param[in] value         the value to assign to all elements of the new DenseVector.
     * @param[in] context   specifies optionally the context where dense vector should reside
     */
    DenseVector( dmemo::DistributionPtr distribution, const ValueType value, hmemo::ContextPtr context = hmemo::ContextPtr() );

    /**
     * @brief creates a replicated DenseVector of the passed size initilized a sequence of values
     *        starting wiht startValue, increased by inc, e.g. [5, 15, 25, 35] with value 5, inc 10
     * 
     * @param[in] size       the size of the new DenseVector.
     * @param[in] startValue the first value of the new DenseVector
     * @param[in] inc        the increment for the sequence of values 
     * @param[in] context    specifies optionally the context where dense vector should reside
     */
    DenseVector( const IndexType size, const ValueType startValue, const ValueType inc, hmemo::ContextPtr context = hmemo::ContextPtr() );

    /**
     * @brief creates a distributed DenseVector of the passed size initilized a sequence of values
     *        starting wiht startValue, increased by inc, e.g. [5, 15, 25, 35] with value 5, inc 10
     * 
     * @param[in] distribution  the distribution to use for the new vector.
     * @param[in] startValue    the first value of the new DenseVector
     * @param[in] inc           the increment for the sequence of values 
     * @param[in] context    specifies optionally the context where dense vector should reside
     */
    DenseVector( dmemo::DistributionPtr distribution, const ValueType startValue, const ValueType inc, hmemo::ContextPtr context = hmemo::ContextPtr() );

    /** Constructor of a replicated vector by replicated C++ array. */

    /**
     * @brief creates a new replicated DenseVector initialized with the passed values.
     *
     * @param[in] size      the size of the new DenseVector.
     * @param[in] values    the values to initialize the new DenseVector with.
     * @param[in] context   specifies optionally the context where dense vector should reside
     */
    template<typename OtherValueType>
    DenseVector( const IndexType size, const OtherValueType* values, hmemo::ContextPtr context = hmemo::ContextPtr() );

    /**
     * Override the default copy constructor to guarantee a deep copy.
     *
     * @param[in] other the dense vector that will be copied
     *
     * The new constructed vector has the same distribution as the input vector.
     */
    DenseVector( const DenseVector<ValueType>& other );

    /**
     * More general constructor that creates a deep copy of an arbitrary vector.
     */

    DenseVector( const Vector& other );

    /**
     * @brief creates a redistributed copy of the passed vector
     *
     * @param[in] other         the vector to take a copy from
     * @param[in] distribution  the distribution to use for the new vector.
     *
     * Must be valid: other.size() == distribution.getGlobalSize()
     */
    DenseVector( const Vector& other, dmemo::DistributionPtr distribution );

    /**
     * @brief creates a distributed DenseVector with given local values.
     *
     * @param[in] localValues   the local values to initialize the new DenseVector with.
     * @param[in] distribution  the distribution the
     */
    DenseVector( const hmemo::_HArray& localValues, dmemo::DistributionPtr distribution );

    /**
     * @brief This constructor creates a vector with the size and values stored
     *        in the file with the given filename.
     *
     * @param[in] filename  the name of the file to read the Vector from.
     *
     * The distribution of the vector is CYCLIC(n) where n is the size of
     * the vector. So only the first processor will hold all values.
     *
     * Note: Only the first processor will read the matrix file.
     */
    DenseVector( const std::string& filename );

    /**
     * @brief creates a DenseVector with the Expression alpha * x.
     *
     * @param[in] expression    alpha * x
     */
    DenseVector( const Expression_SV& expression );

    /**
     * @brief creates a DenseVector with the Expression alpha + x.
     *
     * @param[in] expression    alpha * x + beta
     */
    DenseVector( const Expression_SV_S& expression );

    /**
     *  @brief creates a DenseVector with the Expression alpha * x * Y.
     *
     * @param[in] expression    x * y
     */

    DenseVector( const Expression_VV& expression );


    /**
     *  @brief creates a DenseVector with the Expression alpha * x * Y.
     *
     * @param[in] expression    alpha * x * y
     */

    DenseVector( const Expression_SVV& expression );

    /**
     * @brief creates a DenseVector with the Expression alpha * x + beta * y.
     *
     * @param[in] expression  is alpha * x + beta * y
     */
    DenseVector( const Expression_SV_SV& expression );

    /* --------------------------------------------------------------------- */

    /**
     * @brief creates a DenseVector with the Expression alpha * A * x + beta * y.
     *
     * @param[in] expression     alpha * A * x + beta * y
     */
    DenseVector( const Expression_SMV_SV& expression );

    /**
     * @brief creates a DenseVector with the Expression alpha * x * A + beta * y.
     *
     * @param[in] expression     alpha * x * A + beta * y
     */
    DenseVector( const Expression_SVM_SV& expression );

    /**
     * @brief creates a DenseVector with the Expression alpha * A * x.
     *
     * @param[in] expression     alpha * A * x
     */
    DenseVector( const Expression_SMV& expression );

    /**
     * @brief creates a DenseVector with the Expression alpha * x * A.
     *
     * @param[in] expression     alpha * x * A
     */
    DenseVector( const Expression_SVM& expression );

    /**
     * @brief creates a DenseVector with the Expression A * x.
     *
     * @param[in] expression     A * x
     */
    DenseVector( const Expression_MV& expression );

    /**
     * @brief creates a DenseVector with the Expression x * A.
     *
     * @param[in] expression     x * A
     */
    DenseVector( const Expression_VM& expression );

    /**
     * @brief releases all allocated resources.
     */
    virtual ~DenseVector();

    /** Implementation of pure method Vector::getVectorKind() */

    inline virtual VectorKind getVectorKind() const;

    /** Implememenation of pure routine Vector::allocate. */

    virtual void allocate( dmemo::DistributionPtr distribution );

    /** Implememenation of pure routine Vector::allocate. */

    virtual void allocate( const IndexType n );

    /** Override the default assignment operator.
     *
     *  Note: all other assignment operators are inherited from class Vector.
     */

    DenseVector& operator=( const DenseVector<ValueType>& other );

    /** Reimplement Vector::operator= as otherwise a constructor of DenseVector might be called */

    DenseVector& operator=( const Scalar );

    /**
     * This method initializes a distributed vector with random numbers. 
     * 
     * @param[in] distribution specifies the distribution of the vector
     * @param[in] fillRate for the number of non-zeros
     */
    virtual void setRandom( dmemo::DistributionPtr distribution, const float fillRate = 1.0 );

    /** Implementation of pure method Vector::setSequence */

    virtual void setSequence( const Scalar startValue, const Scalar inc, const IndexType n );

    /** Implementation of pure method Vector::setSequence */

    virtual void setSequence( const Scalar startValue, const Scalar inc, dmemo::DistributionPtr distribution );

    virtual common::scalar::ScalarType getValueType() const;

    /**
     * Implementation of pure method.
     */
    virtual void buildValues( hmemo::_HArray& values ) const;

    /**
     * Implementation of pure method.
     */
    virtual void setValues( const hmemo::_HArray& values );

    /**
     * Implementation of Vector::copy with covariant return type.
     */
    virtual DenseVector* copy() const;

    /**
     * Implementation of Vector::newVector with covariant return type.
     */
    virtual DenseVector* newVector() const;

    //TODO: We either need a none const getLocalValues()
    // or an operator[] with local sematics or both
    // i guess both is the most intuitive because a user that does not request
    // a parallel environment expects the operator[] to exists
    // and for users of a distributed DenseVector the usage of
    // getLocalValues and getHaloValues is more explicite and there for
    // better understandable and less errorprone.
    // Maybe an access proxy would be a nice solution, because with a proxy we
    // can avoid to change the size and other attributes of the HArray
    // mLocalValues.
    /**
     * @brief get a non constant reference to local values of this Dense Vector.
     *
     * @return  a non constant reference to the local values of this.
     */

    utilskernel::LArray<ValueType>& getLocalValues()
    {
        return mLocalValues;
    }

    /**
     * @brief get a constant reference to local values of this Dense Vector.
     *
     * @return  a constant reference to the local values of this.
     */
    const utilskernel::LArray<ValueType>& getLocalValues() const
    {
        return mLocalValues;
    }

    /**
     * @brief Get a reference to the halo temp array of this Dense Vector.
     *
     * @return  a reference to the halo values of this.
     *
     * The halo array can be used as temporary array to keep values of neighbored
     * processors. It avoids reallocation of memory for the values.
     */

    utilskernel::LArray<ValueType>& getHaloValues() const
    {
        return mHaloValues;
    }

    virtual Scalar getValue( IndexType globalIndex ) const;

    void setValue( const IndexType globalIndex, const Scalar value );

    virtual Scalar min() const;

    virtual Scalar max() const;

    virtual Scalar sum() const;

    virtual Scalar l1Norm() const;

    virtual Scalar l2Norm() const;

    virtual Scalar maxNorm() const;

    virtual void swap( Vector& other );

    virtual void writeAt( std::ostream& stream ) const;

    virtual void assign( const Expression_SV_SV& expression );

    virtual void assign( const Expression_SVV& expression );

    virtual void assign( const Expression_SV_S& expression );

    /** Assign this vector with a scalar values, does not change size, distribution. */

    virtual void assign( const Scalar value );

    virtual void add( const Scalar value );

    /** Assign this vector with another vector, inherits size and distribution. */

    virtual void assign( const Vector& other );

    virtual void assign( const hmemo::_HArray& localValues, dmemo::DistributionPtr dist );

    virtual void assign( const hmemo::_HArray& globalValues );

    virtual void buildLocalValues( hmemo::_HArray& localValues ) const;

    virtual Scalar dotProduct( const Vector& other ) const;

    virtual DenseVector& scale( const Vector& other );

    using Vector::prefetch; // prefetch() with no arguments

    virtual void prefetch( const hmemo::ContextPtr location ) const;

    virtual void wait() const;

    virtual void invert();

    virtual void conj();

    virtual void exp();

    virtual void log();

    virtual void floor();

    virtual void ceil();

    virtual void sqrt();

    virtual void sin();

    virtual void cos();

    virtual void tan();

    virtual void atan();

    virtual void powBase( const Vector& other );

    virtual void powExp( const Vector& other );

    virtual void powBase( const Scalar base );

    virtual void powExp( const Scalar exp );

    virtual size_t getMemoryUsage() const;

    virtual void redistribute( dmemo::DistributionPtr distribution );

protected:

    using Vector::mContext;

    SCAI_LOG_DECL_STATIC_LOGGER( logger )

private:

    utilskernel::LArray<ValueType> mLocalValues; //!< my local values of vector

    /** array that might be used to keep halo values of vector, avoids reallocation of memory for halo values */

    mutable utilskernel::LArray<ValueType> mHaloValues; 

public:

    // static methods, variables to register create routine in Vector factory of base class.

    static Vector* create();

    // key for factory

    static VectorCreateKeyType createValue();

    virtual VectorCreateKeyType getCreateValue() const;
};

/* ------------------------------------------------------------------------- */

template<typename ValueType>
Vector::VectorKind DenseVector<ValueType>::getVectorKind() const
{
    return Vector::DENSE;
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
template<typename OtherValueType>
DenseVector<ValueType>::DenseVector( const IndexType size, const OtherValueType* values, hmemo::ContextPtr context )
    : Vector( size, context )
{
    // use LAMA array reference to avoid copy of the raw data
    hmemo::HArrayRef<OtherValueType> valuesArrayRef( size, values );
    // use mContext instead of context to avoid NULL pointer
    utilskernel::HArrayUtils::assign( mLocalValues, valuesArrayRef, mContext );
    // Halo is not used yet
}

} /* end namespace lama */

} /* end namespace scai */
