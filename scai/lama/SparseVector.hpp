/**
 * @file SparseVector.hpp
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
#include <scai/dmemo/Distribution.hpp>
#include <scai/hmemo.hpp>

#include <scai/tasking/SyncToken.hpp>

#include <scai/common/macros/throw.hpp>

// std
#include <fstream>

namespace scai
{

namespace lama
{

// forward declaration required as we do not include DenseVector.hpp here

template<typename ValueType> class DenseVector;

/**
 * @brief SparseVector represents a distributed 1D Vector with elements of type ValueType.
 *
 * In contrary to a dense vector it contains only the non-zero entries.
 *
 * @tparam ValueType the value type for the vector values.
 */
template<typename ValueType>
class COMMON_DLL_IMPORTEXPORT SparseVector:

    public Vector<ValueType>,
    public _Vector::Register<SparseVector<ValueType> >    // register at factory

{
public:

    using _Vector::logger;
    using _Vector::getDistribution;
    using _Vector::getContextPtr;
    using _Vector::readFromFile;
    using _Vector::setDistributionPtr;
    using _Vector::getDistributionPtr;
    using _Vector::size;
    using _Vector::assign;
    using _Vector::prefetch; // prefetch() with no arguments

    using Vector<ValueType>::operator=;
    using Vector<ValueType>::getValueType;
    using Vector<ValueType>::binaryOp;

    /** Default constructor, creates zero-sized vector
     *
     *  @param[in] context  optional, the context to use for this vector.
     */
    SparseVector( hmemo::ContextPtr context = hmemo::Context::getContextPtr() );

    /** Construct a sparse vector of a certain size, all elements are zero */

    SparseVector( const IndexType n, const ValueType zero = 0, hmemo::ContextPtr context = hmemo::Context::getContextPtr() );

    /**
     * @brief creates a distributed SparseVector of the passed global size.
     *
     * @param[in] distribution  the distribution to use for the new vector.
     * @param[in] zero          is the value for all elements
     * @param[in] context       the context to use for the new vector.
     */
    SparseVector( dmemo::DistributionPtr distribution, const ValueType zero = 0, hmemo::ContextPtr context = hmemo::Context::getContextPtr() );

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
        hmemo::ContextPtr context = hmemo::Context::getContextPtr() );

    /**
     * Override the default copy constructor to guarantee a deep copy.
     *
     * @param[in] other the sparse vector that will be copied
     *
     * The new constructed vector has the same distribution as the input vector.
     */
    SparseVector( const SparseVector<ValueType>& other );

    /**
     * @brief Override the default move constructor with appropriate version.
     *
     * @param[in] other the sparse vector from which data is moved
     *
     * The new created object gets its resources from the passed sparse vector.
     */
    SparseVector( SparseVector<ValueType>&& other ) noexcept;

    /** 
     * @brief Construct a vector by a dense vector and zero element.
     *
     * @param[in] other       dense vector to take a copy from
     * @param[in] zeroValue   value that becomes the zero element of the sparse vector.
     */
    explicit SparseVector( const DenseVector<ValueType>& other, const ValueType zeroValue );

    /**
     * @brief creates a distributed SparseVector with given local values.
     *
     * @param[in] distribution  the distribution of the vector
     * @param[in] indexes       local indexes this processor has nonz-zeros values for
     * @param[in] values        values belonging to the corresponding entries
     * @param[in] zero          is the zero value
     * @param[in] context       context to be used for the created sparse vector
     *
     * While indexes and values might have individual values on each processor, the zero value
     * must be the same on all processors
     */
    SparseVector( 
        dmemo::DistributionPtr distribution,
        hmemo::HArray<IndexType> indexes, 
        hmemo::HArray<ValueType> values, 
        const ValueType zero = ValueType( 0 ),
        hmemo::ContextPtr context = hmemo::Context::getContextPtr() );

    /** @brief Short form of SparseVector( dist, indexes, values,  ) if dist is replicated */

    SparseVector( 
        const IndexType n,
        hmemo::HArray<IndexType> indexes, 
        hmemo::HArray<ValueType> values, 
        const ValueType zero = ValueType( 0 ),
        hmemo::ContextPtr context = hmemo::Context::getContextPtr() );

    /**
     * @brief releases all allocated resources.
     */
    virtual ~SparseVector();

    /* ==================================================================== */
    /*   split up member variables                                          */
    /* ==================================================================== */

    void splitUp( dmemo::DistributionPtr& dist,
                  hmemo::HArray<IndexType>& nonZeroIndexes,
                  hmemo::HArray<ValueType>& nonZeroValues,
                  ValueType& zero );

    /** Implementation of pure method _Vector::getVectorKind() */

    inline virtual VectorKind getVectorKind() const;

    /**
     * @brief Implementation of pure method _Vector::isConsistent 
     */
    virtual bool isConsistent() const;

    /** Implememenation of pure routine Vector::allocate. */

    virtual void allocate( dmemo::DistributionPtr distribution );

    /** Implememenation of pure routine _Vector::allocate. */

    virtual void allocate( const IndexType n );

    /** Query the zero value, i.e. default value at positions not in nonZeroIndexes. */
    
    inline ValueType getZero() const;

    /** Override the default assignment operator.  */

    SparseVector& operator=( const SparseVector<ValueType>& other );

    /** Override the default move assignment operator. */

    SparseVector& operator=( SparseVector<ValueType>&& other ) noexcept;

    /** Implementation of pure method _Vector::assign 
     *
     *  uses metaprogramming to call assignImpl with actual type and kind
     */
    virtual void assign( const _Vector& other );

    template<typename OtherValueType>
    void assignImpl( const Vector<OtherValueType>& other );

    template<typename OtherValueType>
    void assignSparse( const SparseVector<OtherValueType>& other );

    template<typename OtherValueType>
    void assignDense( const DenseVector<OtherValueType>& other );

    /**
     * Implementation of _Vector::fillRandom for sparse vectors.
     *
     * This method is only available to keep consistency with the dense vector.
     * If applied to a sparse vector it would not remain really sparse. 
     * Please use  setSparseRandom to generate a sparse random vector.
     */
    virtual void fillRandom( const IndexType bound );

    /**
     *  This method generates a random sparse vector with a certain fillRate.
     * 
     *  @param[in] bound     draw random non-zero values in range of $[0,bound]$.
     *  @param[in] fillRate  probablity for a non-zero entry
     */
    virtual void fillSparseRandom( const float fillRate, const IndexType bound );

    /**
     * @brief Implementation of pure method _Vector::buildLocalValues.
     * 
     * For a sparse vector this routine builds a dense local part (COPY) or
     * it scatters its values in the localValues array corresponding to op.
     */
    virtual void buildLocalValues(
        hmemo::_HArray& localValues,
        const common::BinaryOp op = common::BinaryOp::COPY,
        hmemo::ContextPtr prefLoc = hmemo::ContextPtr() ) const;

    /**
     * @brief Implementation of pure method _Vector::gatherLocalValues.
     * 
     * For a sparse vector each of the required indexes must be searched in the
     * nonZeroIndexes. The corresponding position can be used to get the required value.
     */
    virtual void gatherLocalValues(
        hmemo::_HArray& localValues,
        const hmemo::HArray<IndexType>& localIndexes,
        const common::BinaryOp op = common::BinaryOp::COPY,
        hmemo::ContextPtr prefLoc = hmemo::ContextPtr() ) const;

    /** Get the array with the non-zero values. */

    const hmemo::HArray<ValueType>& getNonZeroValues() const;

    /** Get the array with the non-zero indexes. */

    const hmemo::HArray<IndexType>& getNonZeroIndexes() const;

    /**
     * Implementation of pure method _Vector::setDenseValues.
     */
    virtual void setDenseValues( const hmemo::_HArray& values );

    /**
     *  Typed version of setDenseValues, i.e. _HArray is now HArray<OtherValueType>
     */
    template<typename OtherValueType>
    void setDenseValuesImpl( const hmemo::HArray<OtherValueType>& values );

    /**
     * Implementation of pure method _Vector::fillSparseData
     */
    virtual void fillSparseData( 
        const hmemo::HArray<IndexType>& nonZeroIndexes, 
        const hmemo::_HArray& nonZeroValues,
        const common::BinaryOp op );

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
     * Implementation of _Vector::copy with covariant return type.
     */
    virtual SparseVector* copy() const;

    /**
     * Implementation of _Vector::newVector with covariant return type.
     */
    virtual SparseVector* newVector() const;

    /** @brief Implementation of pure method Vector<ValueType>::getValue */

    virtual ValueType getValue( IndexType globalIndex ) const;

    /** @brief Implementation of pure method Vector<ValueType>::setValue */

    virtual void setValue( const IndexType globalIndex, const ValueType value );

    virtual ValueType min() const;

    virtual ValueType max() const;

    virtual ValueType sum() const;

    virtual RealType<ValueType> l1Norm() const;

    virtual RealType<ValueType> l2Norm() const;

    /** Implementation of pure method _Vector::maxNorm */

    virtual RealType<ValueType> maxNorm() const;

    /** Implementation of pure method Vector<ValueType>::maxDiffNorm */

    virtual RealType<ValueType> maxDiffNorm( const Vector<ValueType>& other ) const;

    /** Implementation of pure method Vector<ValueType>::all */

    virtual bool all( common::CompareOp op, const ValueType alpha ) const;

    /** Implementation of pure method Vector<ValueType>::all */

    virtual bool all( common::CompareOp op, const Vector<ValueType>& x ) const;

    virtual void swap( _Vector& other );

    virtual void writeAt( std::ostream& stream ) const;

    /** Implementation for pure method Vector<ValueType>::setScalar for sparse vector 
     *
     *  Note: value becomes the zero element of the sparse vector.
     */
    virtual void setScalar( const ValueType& alpha );

    /** Implementation of pure method Vector<ValueType>::unaryOp for sparse vector. */

    void unaryOp( const Vector<ValueType>& x, const common::UnaryOp op );

    /** Implementation of pure method Vector<ValueType>::binaryOp for sparse vector. */

    void binaryOp( const Vector<ValueType>& x, const common::BinaryOp op, const Vector<ValueType>& y );

    /** Implementation of pure method Vector<ValueType>::binaryOpScalar for sparse vector. */

    void binaryOpScalar( const Vector<ValueType>& x, const ValueType& alpha, const common::BinaryOp op, const bool swap );

    /** Implementation of pure method Vector<ValueType>::selectComplexPart for sparse vector. */

    virtual void selectComplexPart( Vector<RealType<ValueType> >& x, const common::ComplexPart part ) const;

    /** Implementation of pure method Vector<ValueType>::buildComplex for sparse vector. */

    virtual void buildComplex( const Vector<RealType<ValueType> >& x, const Vector<RealType<ValueType> >& y );

    /** Implementation of pure method Vector<ValueType>::vectorPlusVector */

    virtual void vectorPlusVector( const ValueType& alpha, const Vector<ValueType>& x, const ValueType& beta, const Vector<ValueType>& y );

    /** Implmentation of vectorPlusVector for sparse vectors of same type */

    void vectorPlusVectorImpl( const ValueType alpha, const SparseVector<ValueType>& x, 
                               const ValueType beta, const SparseVector<ValueType>& y );

    /** Implementation of pure method Vector<ValueType>::vectorTimesVector */

    virtual void vectorTimesVector( const ValueType& alpha, const Vector<ValueType>& x, const Vector<ValueType>& y );

    /** Implementation of pure method Vector<ValueType>::vectorPlusScalar */

    virtual void vectorPlusScalar( const ValueType& alpha, const Vector<ValueType>& x, const ValueType& beta );

    /** Implementation of pure method Vector::dotProduct */

    virtual ValueType dotProduct( const Vector<ValueType>& other ) const;

    /** Implementation of pure method _Vector::setVector */

    virtual void setVector( const _Vector& other, const common::BinaryOp op, const bool swapArgs = false );

    virtual void prefetch( const hmemo::ContextPtr location ) const;

    virtual void wait() const;

    virtual size_t getMemoryUsage() const;

    /** Implementation of pure method _Vector::redistribute */

    virtual void redistribute( dmemo::DistributionPtr distribution );

    /** Implementation of pure method _Vector::redistribute */

    virtual void redistribute( const dmemo::RedistributePlan& redistributor );

    /** Implementation of pure method _Vector::resize */

    virtual void resize( const dmemo::DistributionPtr distribution );

private:

    /** Help predicate, checks for same sparse pattern on local part only, no global comm. */

    bool hasSamePattern( const SparseVector<ValueType>& other ) const;

    ValueType dotProductLocalAny( const Vector<ValueType>& other ) const;

    ValueType dotProductLocalSparse( const SparseVector<ValueType>& other ) const;

    ValueType dotProductLocalDense( const DenseVector<ValueType>& other ) const;

    /** Help routine that mults a sparse vector ( zero element is 0 ) with a dense vector */

    void binaryOpMult0D( const SparseVector<ValueType>& x, const DenseVector<ValueType>& y );

    /** Help routine that mults a sparse vector ( zero element is 0 ) with a any other vector */

    void binaryOpMult0( const SparseVector<ValueType>& x, const Vector<ValueType>& y );

    void binaryOpAny( const Vector<ValueType>& x, const common::BinaryOp op, const Vector<ValueType>& y );

    void binaryOpSparse( const SparseVector<ValueType>& x, const common::BinaryOp op, const SparseVector<ValueType>& y );

    /** Help routine for binary operation of two sparse vectors with same pattern */

    void binaryOpSparseSamePattern( const SparseVector<ValueType>& x, const common::BinaryOp op, const SparseVector<ValueType>& y );

    /** Help routine for binary operation of two sparse vectors with different pattern */

    void binaryOpSparseNewPattern( const SparseVector<ValueType>& x, const common::BinaryOp op, const SparseVector<ValueType>& y );

    hmemo::HArray<IndexType> mNonZeroIndexes;  //!< my local indexes for non-zero values
    hmemo::HArray<ValueType> mNonZeroValues;   //!< my local non-zero values

    ValueType mZeroValue;    //!< ZERO element of vector, can be explicity set, defaults to 0

    /** array that might be used to keep halo values of vector, avoids reallocation of memory for halo values */

    mutable hmemo::HArray<ValueType> mHaloValues;

    /** Implementation of _Vector::writeLocalToFile */

    virtual void writeLocalToFile( FileIO& file ) const;

    /** Implementation of _Vector::readFromFile */

    virtual void readFromFile( FileIO& file );

    /** Implementation of _Vector::clearValues */

    virtual void clearValues();

public:

    // static methods, variables to register create routine in Vector factory of base class.

    static _Vector* create();

    // key for factory

    static VectorCreateKeyType createValue();
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

    Vector<ValueType>( size, context )

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
    setScalar( zeroValue );

    // use LAMA array reference to avoid copy of the raw data

    hmemo::HArrayRef<OtherValueType> values( nnz, nonZeroValues );
    hmemo::HArrayRef<IndexType> indexes( nnz, nonZeroIndexes );

    // values at same index will be replaced

    fillSparseData( indexes, values, common::BinaryOp::COPY );
}

template<typename ValueType>
VectorKind SparseVector<ValueType>::getVectorKind() const
{
    return VectorKind::SPARSE;
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
ValueType SparseVector<ValueType>::getZero() const
{
    return mZeroValue;
}

/** 
 *  @brief Free function that returns a sparse vector of a given size initialized with zero value.
 * 
 *  @tparam    ValueType  is the component type of the vector
 *  @param[in] n          specifies the size of the vector                             
 *  @param[in] zero       is the zero element of the sparse vector
 *  @param[in] ctx        Context that is used for the filling and the generated vector
 *  @returns              a new sparse vector with the specified size
 *
 *  \code
 *     const auto v = sparseVector<double>( n, 1.0 );
 *  \endcode
 */
template<typename ValueType>
SparseVector<ValueType> sparseVector(
    const IndexType n,
    ValueType zero = 0, 
    hmemo::ContextPtr ctx = hmemo::Context::getContextPtr() )
{
    return SparseVector<ValueType>( n, zero, ctx );
}

/** 
 *  @brief Free function that returns a sparse vector with a given distributon initialized with zero value.
 * 
 *  @tparam    ValueType    is the component type of the vector
 *  @param[in] distribution specifies size and mapping of the vector
 *  @param[in] zero         is the zero element of the sparse vector
 *  @param[in] ctx          Context that is used for the filling and the generated vector
 *  @returns                a new sparse vector with the specified distributiondistribution.
 *
 *  \code
 *     const auto v = sparseVector<double>( blockDistribution( n ), 1.0 );
 *  \endcode
 */
template<typename ValueType>
SparseVector<ValueType> sparseVector(
    dmemo::DistributionPtr distribution,
    ValueType zero = 0,
    hmemo::ContextPtr ctx = hmemo::Context::getContextPtr() )
{
    return SparseVector<ValueType>( distribution, zero, ctx );
}

} /* end namespace lama */

} /* end namespace scai */
