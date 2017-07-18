/**
* @file SparseVector.hpp
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
* @brief Definition of template class that stands for a sparse vector of a certain type.
* @author Thomas Brandes
* @date 16.01.2017
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
* @brief Common base class for all sparse vectors to deal with untyped sparse vectors.
*
* A sparse vector is a distributed one-dimensional array where only non-zero values are explicitly stored.
* This base class provides common methods for typed sparse vectors and allows access to the
* data via untyped heterogeneous arrays.
*/

class COMMON_DLL_IMPORTEXPORT _SparseVector :

public Vector

{

public:

/** Implementation of pure method Vector::getVectorKind() */

inline virtual VectorKind getVectorKind() const;

/**
* @brief get a constant reference to local non-zero values of this Sparse Vector.
*
* @return  a constant reference to the local non-zero values of this dense vector.
*/
virtual const hmemo::_HArray& getNonZeroValues() const = 0;

virtual const hmemo::HArray<IndexType>& getNonZeroIndexes() const = 0;

/** Query the zero value, i.e. default value at positions not in nonZeroIndexes. */

virtual Scalar getZero() const = 0;

/**
* @brief Implementation of pure method Vector::isConsistent 
*/
virtual bool isConsistent() const;

/**
     * @brief Create a new sparse vector of a certain type 
     *
     * @param type is the value type of the vector
     * @return new allocated SparseVector of the given type 
     * @throw  common::Exception if no dense vector of this type is registered in factory.
     */

    static _SparseVector* create( common::scalar::ScalarType type );

    // make operators and methods of Vector visible for _SparseVector

    // make operators and methods of Vector visible for _SparseVector

    using Vector::operator=;

protected:

    // All constructors here just call the corresponing constructors of Vector 

    _SparseVector( const IndexType n );

    _SparseVector( const IndexType n, hmemo::ContextPtr context );

    _SparseVector( const dmemo::DistributionPtr dist );

    _SparseVector( const dmemo::DistributionPtr dist, hmemo::ContextPtr context );

    _SparseVector( const _SparseVector& other );

    _SparseVector( const Vector& other );

    /** Common logger for all types of sparse vectors. */

    SCAI_LOG_DECL_STATIC_LOGGER( logger )
};

/**
 * @brief SparseVector represents a distributed 1D Vector with elements of type ValueType.
 *
 * In contrary to a dense vector it contains only the non-zero entries.
 *
 * @tparam ValueType the value type for the vector values.
 */
template<typename ValueType>
class COMMON_DLL_IMPORTEXPORT SparseVector:

    public _SparseVector,

    public Vector::Register<SparseVector<ValueType> >    // register at factory

{
public:

    /** Default constructor, creates zero-sized replicated vector */

    SparseVector();

    /** Construct a sparse vector of a certain size, all elements are zero */

    explicit SparseVector( const IndexType n );

    /** Construct a sparse vector of a certain size, all elements are zero, with a certain context */

    SparseVector( const IndexType n, hmemo::ContextPtr context );

    /**
     * @brief creates a not initialized distributed SparseVector of the passed global size.
     *
     * @param[in] context  the context to use for the new vector.
     */

    explicit SparseVector( hmemo::ContextPtr context );

    /**
     * @brief creates a not initialized distributed SparseVector of the passed global size.
     *
     * @param[in] distribution  the distribution to use for the new vector.
     */
    explicit SparseVector( dmemo::DistributionPtr distribution );

    /**
     * @brief creates a not initialized distributed SparseVector of the passed global size.
     *
     * @param[in] distribution  the distribution to use for the new vector.
     * @param[in] context  the context to use for the new vector.
     */
    SparseVector( dmemo::DistributionPtr distribution, hmemo::ContextPtr context );

    /**
     * @brief creates a replicated SparseVector of the passed size with value as its ZERO element
     *
     * @param[in] size  the size of the new SparseVector.
     * @param[in] value the value becomes ZERO element of the sparse vector, i.e. all elements have this values
     * @param[in] context specifies optionally the context where dense vector should reside
     */
    SparseVector( const IndexType size, const ValueType value, hmemo::ContextPtr context = hmemo::ContextPtr() );

    /**
     * @brief Constrcutor of sparse vector with raw arrays.
     */
    template<typename OtherValueType>
    SparseVector(
        const IndexType size,
        const IndexType nnz,
        const IndexType nonZeroIndexes[],
        const OtherValueType nonZeroValues[],
        const OtherValueType zeroValue = OtherValueType( 0 ),
        hmemo::ContextPtr context = hmemo::ContextPtr() );

    /**
     * Override the default copy constructor to guarantee a deep copy.
     *
     * @param[in] other the sparse vector that will be copied
     *
     * The new constructed vector has the same distribution as the input vector.
     */
    SparseVector( const SparseVector<ValueType>& other );

    /**
     * More general constructor that creates a deep copy of an arbitrary vector.
     *
     * The explicit specifier avoids implict conversions as the following example shows.
     *
     * \code
     *     subroutine sub( const SparseVector<float>& v );
     *     ...
     *     SparseVector<float> vf;
     *     SparseVector<double> vd;
     *     sub( vf );               // that is okay
     *     sub( vd );               // compile error to avoid implicit conversions
     * \endcode
     */
    explicit SparseVector( const Vector& other );

    /**
     * @brief creates a redistributed copy of the passed vector
     *
     * @param[in] other         the vector to take a copy from
     * @param[in] distribution  the distribution to use for the new vector.
     *
     * Must be valid: other.size() == distribution.getGlobalSize()
     */
    explicit SparseVector( const Vector& other, dmemo::DistributionPtr distribution );

    /**
     * @brief creates a distributed SparseVector with given local values.
     *
     * @param[in] distribution  the distribution of the vector
     * @param[in] indexes       local indexes this processor has nonz-zeros values for
     * @param[in] values        values belonging to the corresponding entries
     * @param[in] zero          is the zero value
     *
     * While indexes and values might have individual values on each processor, the zero value
     * must be exactly the same on all processors
     */
    SparseVector( 
        dmemo::DistributionPtr distribution,
        const hmemo::HArray<IndexType>& indexes, 
        const hmemo::_HArray& values, 
        const Scalar zero = Scalar( 0 ) );

    /**
     * @brief creates a distributed SparseVector with given local values.
     *
     * @param[in] distribution  the distribution of the vector
     * @param[in] localValues   local values
     *
     * This constructor is only available to be more conform with DenseVector.
     * The zero value is assumed to be 0. Only non-zero entries are stored explicitly.
     */
    SparseVector( dmemo::DistributionPtr distribution, const hmemo::_HArray& localValues );

    /** @brief Short form of SparseVector( dist, localValues ) if dist is replicated */

    explicit SparseVector( const hmemo::_HArray& localValues );

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
    explicit SparseVector( const std::string& filename );

    /**
     * @brief creates a SparseVector with the Expression alpha * x.
     *
     * @param[in] expression    alpha * x
     */
    explicit SparseVector( const Expression_SV& expression );

    /**
     * @brief creates a SparseVector with the Expression alpha + x.
     *
     * @param[in] expression    alpha * x + beta
     */
    explicit SparseVector( const Expression_SV_S& expression );

    /**
     *  @brief creates a SparseVector with the Expression alpha * x * Y.
     *
     * @param[in] expression    x * y
     */

    explicit SparseVector( const Expression_VV& expression );


    /**
     *  @brief creates a SparseVector with the Expression alpha * x * Y.
     *
     * @param[in] expression    alpha * x * y
     */

    explicit SparseVector( const Expression_SVV& expression );

    /**
     * @brief creates a SparseVector with the Expression alpha * x + beta * y.
     *
     * @param[in] expression  is alpha * x + beta * y
     */
    explicit SparseVector( const Expression_SV_SV& expression );

    /* --------------------------------------------------------------------- */

    /**
     * @brief creates a SparseVector with the Expression alpha * A * x + beta * y.
     *
     * @param[in] expression     alpha * A * x + beta * y
     */
    explicit SparseVector( const Expression_SMV_SV& expression );

    /**
     * @brief creates a SparseVector with the Expression alpha * x * A + beta * y.
     *
     * @param[in] expression     alpha * x * A + beta * y
     */
    explicit SparseVector( const Expression_SVM_SV& expression );

    /**
     * @brief creates a SparseVector with the Expression alpha * A * x.
     *
     * @param[in] expression     alpha * A * x
     */
    explicit SparseVector( const Expression_SMV& expression );

    /**
     * @brief creates a SparseVector with the Expression alpha * x * A.
     *
     * @param[in] expression     alpha * x * A
     */
    explicit SparseVector( const Expression_SVM& expression );

    /**
     * @brief creates a SparseVector with the Expression A * x.
     *
     * @param[in] expression     A * x
     */
    explicit SparseVector( const Expression_MV& expression );

    /**
     * @brief creates a SparseVector with the Expression x * A.
     *
     * @param[in] expression     x * A
     */
    explicit SparseVector( const Expression_VM& expression );

    /**
     * @brief releases all allocated resources.
     */
    virtual ~SparseVector();

    /** Implememenation of pure routine Vector::allocate. */

    virtual void allocate( dmemo::DistributionPtr distribution );

    /** Implememenation of pure routine Vector::allocate. */

    virtual void allocate( const IndexType n );

    /** Implementation of pure method _SparseVector::getZero */

    inline Scalar getZero() const;

    /** Override the default assignment operator.  */

    SparseVector& operator=( const SparseVector<ValueType>& other );

    // All other assignment operators are inherited from class Vector, but using is required

    using Vector::operator=;
    using Vector::assign;

    /**
     * This method initializes a distributed vector with random numbers.
     *
     * @param[in] distribution specifies the distribution of the vector
     * @param[in] fillRate for the number of non-zeros
     */
    virtual void setRandom( dmemo::DistributionPtr distribution, const float fillRate = 1.0 );

    /** Implementation of Vector::getValueType */

    virtual common::scalar::ScalarType getValueType() const;

    /**
     * @brief Implementation of pure method Vector::buildLocalValues.
     * 
     * For a sparse vector this routine builds a dense local part (COPY) or
     * it scatters its values in the localValues array corresponding to op.
     */
    virtual void buildLocalValues(
        hmemo::_HArray& localValues,
        const common::binary::BinaryOp op = common::binary::COPY,
        hmemo::ContextPtr prefLoc = hmemo::ContextPtr() ) const;

    /**
     * @brief Implementation of pure method Vector::gatherLocalValues.
     * 
     * For a sparse vector each of the required indexes must be searched in the
     * nonZeroIndexes. The corresponding position can be used to get the required value.
     */
    virtual void gatherLocalValues(
        hmemo::_HArray& localValues,
        const hmemo::HArray<IndexType>& localIndexes,
        const common::binary::BinaryOp op = common::binary::COPY,
        hmemo::ContextPtr prefLoc = hmemo::ContextPtr() ) const;

    /** Implementation of _SpaseVector::getNonZeroValues */

    virtual const hmemo::HArray<ValueType>& getNonZeroValues() const;

    /** Implementation of _SpaseVector::getNonZeroIndexes */

    virtual const hmemo::HArray<IndexType>& getNonZeroIndexes() const;

    /**
     * Implementation of pure method Vector::setDenseValues.
     */
    virtual void setDenseValues( const hmemo::_HArray& values );

    /**
     * Implementation of pure method Vector::setSparseValues.
     */
    virtual void setSparseValues(
        const hmemo::HArray<IndexType>& nonZeroIndexes,
        const hmemo::_HArray& nonZeroValues,
        const Scalar zeroValue = Scalar( 0 ) );

    /**
     * Setting sparse values with raw data 
     *
     *   @param[in] nnz  number of non-zero values 
     *   @param[in] nonZeroIndexes array with positions of non-zero values
     *   @param[in] nonZeroValues array with the values at the positions
     *   @param[in] zeroValue is the value for all other array elements
     *
     *   Note: the non zero entries will be sorted by position in ascending order.
     */
    template<typename OtherValueType>
    void setSparseValues( 
        const IndexType nnz, 
        const IndexType nonZeroIndexes[], 
        const OtherValueType nonZeroValues[],
        const OtherValueType zeroValue = OtherValueType( 0 ) );

    /**
     * Set local sparse values of the sparse vector by swapping with existing arrays.
     *
     *  @param[in, out] nonZeroIndexes in the new indexes, out the old indexes
     *  @param[in, out] nonZeroValues in the new values, out the old values
     *
     * The new sparse entries will be sorted in any case.
     */
    void swapSparseValues( hmemo::HArray<IndexType>& nonZeroIndexes, hmemo::HArray<ValueType>& nonZeroValues );

    /**
     * Implementation of Vector::copy with covariant return type.
     */
    virtual SparseVector* copy() const;

    /**
     * Implementation of Vector::newVector with covariant return type.
     */
    virtual SparseVector* newVector() const;

    virtual Scalar getValue( IndexType globalIndex ) const;

    virtual void setValue( const IndexType globalIndex, const Scalar value );

    virtual Scalar min() const;

    virtual Scalar max() const;

    virtual Scalar sum() const;

    virtual Scalar l1Norm() const;

    virtual Scalar l2Norm() const;

    /** Implementation of pure method Vector::maxNorm */

    virtual Scalar maxNorm() const;

    /** Implementation of pure method Vector::maxDiffNorm */

    virtual Scalar maxDiffNorm( const Vector& other ) const;

    virtual void swap( Vector& other );

    virtual void writeAt( std::ostream& stream ) const;

    /** Implementation of pure method Vector::vectorPlusVector */

    virtual void vectorPlusVector( const Scalar& alphaS, const Vector& x, const Scalar& betaS, const Vector& y );

    /** Implmentation of vectorPlusVector for sparse vectors of same type */

    void vectorPlusVectorImpl( const ValueType alpha, const SparseVector<ValueType>& x, 
                               const ValueType beta, const SparseVector<ValueType>& y );

    /** Implementation of pure method Vector::vectorTimesVector */

    virtual void vectorTimesVector( const Scalar& alphaS, const Vector& x, const Vector& y );

    /** Implementation of pure method Vector::vectorPlusScalar */

    virtual void vectorPlusScalar( const Scalar& alphaS, const Vector& x, const Scalar& betaS );

    /** Assign this vector with a scalar values, does not change size, distribution. */

    virtual void assign( const Scalar value );

    /** Implementation of pure method Vector::dotProduct */

    virtual Scalar dotProduct( const Vector& other ) const;

    /** Implementation of pure method Vector::setVector */

    virtual void setVector( const Vector& other, const common::binary::BinaryOp op, const bool swapArgs = false );

    /** Implementation of pure method Vector::setScalar */

    virtual void setScalar( const Scalar value, common::binary::BinaryOp op, const bool swapArgs = false );

    /** Implementation of pure method Vector::applyUnary */

    virtual void applyUnary( common::unary::UnaryOp op );

    using Vector::prefetch; // prefetch() with no arguments

    virtual void prefetch( const hmemo::ContextPtr location ) const;

    virtual void wait() const;

    virtual size_t getMemoryUsage() const;

    virtual void redistribute( dmemo::DistributionPtr distribution );

protected:

    using Vector::mContext;

    /** Help routine for binary operation of two sparse vectors */

    void binOpSparse( const _SparseVector& other, const common::binary::BinaryOp op, bool swapArgs );

private:

    utilskernel::LArray<IndexType> mNonZeroIndexes;  //!< my local indexes for non-zero values
    utilskernel::LArray<ValueType> mNonZeroValues;   //!< my local non-zero values

    ValueType mZeroValue;    //!< ZERO element of vector, can be explicity set, defaults to 0

    /** array that might be used to keep halo values of vector, avoids reallocation of memory for halo values */

    mutable utilskernel::LArray<ValueType> mHaloValues;

    /** Implementation of Vector::writeLocalToFile */

    virtual void writeLocalToFile(
        const std::string& fileName,
        const std::string& fileType,
        const common::scalar::ScalarType dataType,
        const FileIO::FileMode fileMode ) const;

    /** Implementation of Vector::readLocalFromFile */

    virtual IndexType readLocalFromFile( const std::string& fileName, const IndexType first = 0, const IndexType size = nIndex );

    /** Implementation of Vector::clearValues */

    virtual void clearValues();

public:

    // static methods, variables to register create routine in Vector factory of base class.

    static Vector* create();

    // key for factory

    static VectorCreateKeyType createValue();

    virtual VectorCreateKeyType getCreateValue() const;
};

/* ------------------------------------------------------------------------- */
/*  Implementation of inline constructors                                    */
/* ------------------------------------------------------------------------- */

template<typename ValueType>
template<typename OtherValueType>
SparseVector<ValueType>::SparseVector(
    const IndexType size,
    const IndexType nnz,
    const IndexType nonZeroIndexes[],
    const OtherValueType nonZeroValues[],
    const OtherValueType zeroValue,
    hmemo::ContextPtr context ) :

    _SparseVector( size, context )

{
    setSparseValues( nnz, nonZeroIndexes, nonZeroValues, zeroValue );
}

/* ------------------------------------------------------------------------- */
/*  Implementation of inline methods                                         */
/* ------------------------------------------------------------------------- */

template<typename ValueType>
template<typename OtherValueType>
void SparseVector<ValueType>::setSparseValues(
    const IndexType nnz,
    const IndexType nonZeroIndexes[],
    const OtherValueType nonZeroValues[],
    const OtherValueType zeroValue )
{
    // use LAMA array reference to avoid copy of the raw data

    hmemo::HArrayRef<OtherValueType> values( nnz, nonZeroValues );
    hmemo::HArrayRef<IndexType> indexes( nnz, nonZeroIndexes );

    setSparseValues( indexes, values, zeroValue );
}

Vector::VectorKind _SparseVector::getVectorKind() const
{
    return Vector::SPARSE;
}

template<typename ValueType>
const hmemo::HArray<ValueType>& SparseVector<ValueType>::getNonZeroValues() const
{
    return mNonZeroValues;
}

template<typename ValueType>
const hmemo::HArray<IndexType>& SparseVector<ValueType>::getNonZeroIndexes() const
{
    return mNonZeroIndexes;
}

template<typename ValueType>
Scalar SparseVector<ValueType>::getZero() const
{
    return Scalar( mZeroValue );
}

} /* end namespace lama */

} /* end namespace scai */
