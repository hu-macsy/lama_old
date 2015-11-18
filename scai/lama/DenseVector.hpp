/**
 * @file DenseVector.hpp
 *
 * @license
 * Copyright (c) 2009-2015
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
 * @brief DenseVector.hpp
 * @author Jiri Kraus
 * @date 22.02.2011
 */

#pragma once

// for dll_import
#include <scai/common/config.hpp>

// base classes
#include <scai/lama/Vector.hpp>

// local library
#include <scai/lama/LAMAArrayUtils.hpp>
#include <scai/lama/distribution/Distribution.hpp>
#include <scai/lama/distribution/Halo.hpp>

#include <scai/lama/io/mmio.hpp>
#include <scai/lama/io/FileType.hpp>
#include <scai/lama/io/XDRFileStream.hpp>

// internal scai libraries
#include <scai/hmemo.hpp>

#include <scai/tasking/SyncToken.hpp>

#include <scai/common/macros/throw.hpp>

// std
#include <fstream>

using namespace scai::tasking;

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

    /** Default constructor, creates replicated 0 vector */

    DenseVector();

    explicit DenseVector( hmemo::ContextPtr context );

    /**
     * @brief creates a not initialized distributed DenseVector of the passed global size.
     *
     * @param[in] distribution  the distribution to use for the new vector.
     */
    explicit DenseVector( DistributionPtr distribution );

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
    DenseVector( DistributionPtr distribution, const ValueType value, hmemo::ContextPtr context = hmemo::ContextPtr() );

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
    DenseVector( const Vector& other, DistributionPtr distribution );

    /**
     * @brief creates a distributed DenseVector with given local values.
     *
     * @param[in] localValues   the local values to initialize the new DenseVector with.
     * @param[in] distribution  the distribution the
     */
    DenseVector( const hmemo::ContextArray& localValues, DistributionPtr distribution );

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
    DenseVector( const Expression<Scalar,Expression<Matrix,Vector,Times>,Times>& expression );

    /**
     * @brief creates a DenseVector with the Expression alpha * x * A.
     *
     * @param[in] expression     alpha * x * A
     */
    DenseVector( const Expression<Scalar,Expression<Vector,Matrix,Times>,Times>& expression );

    /**
     * @brief creates a DenseVector with the Expression A * x.
     *
     * @param[in] expression     A * x
     */
    DenseVector( const Expression<Matrix,Vector,Times>& expression );

    /**
     * @brief creates a DenseVector with the Expression x * A.
     *
     * @param[in] expression     x * A
     */
    DenseVector( const Expression<Vector,Matrix,Times>& expression );

    /**
     * @brief releases all allocated resources.
     */
    virtual ~DenseVector();

    /** Allocate a dense vector with a certain distribution, values are undefined. */

    void allocate( DistributionPtr distribution );

    /** Override the default assignment operator.
     *
     *  Note: all other assignment operators are inherited from class Vector.
     */

    DenseVector& operator=( const DenseVector<ValueType>& other );

    /** Reimplement scalar assignment as otherwise type conversion rules do not apply. */

    DenseVector& operator=( const Scalar );

    /**
     * This method implements Vector::readFromFile.
     *
     * The distribution of the vector is CYCLIC(n) where n is the size of
     * the vector. So only the first processor will hold all values.
     *
     * Note: Only the first processor will read the matrix file.
     */
    void readFromFile( const std::string& filename );

    virtual common::ScalarType getValueType() const;

    /**
     * Implementation of pure method.
     */
    virtual void buildValues( hmemo::ContextArray& values ) const;

    /**
     * Implementation of pure method.
     */
    virtual void setValues( const hmemo::ContextArray& values );

    /**
     * Implementation of Vector::clone with covariant return type.
     */

    virtual DenseVector* clone() const;

    /**
     * Implementation of Vector::clone with covariant return type.
     */
    virtual DenseVector* clone( DistributionPtr distribution ) const;

    /**
     * Implementation of Vector::copy with covariant return type.
     */

    virtual DenseVector* copy() const;

    //TODO: We either need a none const getLocalValues()
    // or an operator[] with local sematics or both
    // i guess both is the most intuitive because a user that does not request
    // a parallel environment expects the operator[] to exists
    // and for users of a distributed DenseVector the usage of
    // getLocalValues and getHaloValues is more explicite and there for
    // better understandable and less errorprone.
    // Maybe an access proxy would be a nice solution, because with a proxy we
    // can avoid to change the size and other attributes of the LAMAArray
    // mLocalValues.
    /**
     * @brief get a non constant reference to local values of this Dense Vector.
     *
     * @return  a non constant reference to the local values of this.
     */

    hmemo::LAMAArray<ValueType>& getLocalValues()
    {
        return mLocalValues;
    }

    /**
     * @brief get a constant reference to local values of this Dense Vector.
     *
     * @return  a constant reference to the local values of this.
     */
    const hmemo::LAMAArray<ValueType>& getLocalValues() const
    {
        return mLocalValues;
    }

    /**
     * @brief get a reference to halo values of this Dense Vector.
     *
     * @return  a reference to the halo values of this.
     *
     * Note: halo of a vector can also be used for writes in case of const vectors.
     */
    hmemo::LAMAArray<ValueType>& getHaloValues() const
    {
        return mHaloValues;
    }

    /**
     * @brief update the halo values according to the passed Halo.
     *
     * @param[in] halo  the halo which describes which remote values should be put into the halo cache.
     */
    void updateHalo( const Halo& halo ) const;

    /**
     * @brief update the halo values according to the passed Halo asynchronously.
     *
     * @param[in] halo  the halo which describes which remote values should be put into the halo cache.
     * @return          a SyncToken which can be used to synchronize to the asynchronous update.
     */
    SyncToken* updateHaloAsync( const Halo& halo ) const;

    virtual Scalar getValue( IndexType globalIndex ) const;

    virtual Scalar min() const;

    virtual Scalar max() const;

    virtual Scalar l1Norm() const;

    virtual Scalar l2Norm() const;

    virtual Scalar maxNorm() const;

    static void vectorPlusVector(
        scai::hmemo::ContextPtr prefContext,
        scai::hmemo::LAMAArray<ValueType>& result,
        const ValueType alpha,
        const scai::hmemo::LAMAArray<ValueType>& x,
        const ValueType beta,
        const scai::hmemo::LAMAArray<ValueType>& y );

    static SyncToken* vectorPlusVectorAsync(
        hmemo::ContextPtr prefContext,
        hmemo::LAMAArray<ValueType>& result,
        const ValueType alpha,
        const hmemo::LAMAArray<ValueType>& x,
        const ValueType beta,
        const hmemo::LAMAArray<ValueType>& y );

    virtual void swap( Vector& other );

    virtual void writeAt( std::ostream& stream ) const;

    virtual void assign( const Expression_SV_SV& expression );

    /** Assign this vector with a scalar values, does not change size, distribution. */

    virtual void assign( const Scalar value );

    /** Assign this vector with another vector, inherits size and distribution. */

    virtual void assign( const Vector& other );

    virtual void assign( const hmemo::ContextArray& localValues, DistributionPtr dist );

    virtual void buildLocalValues( hmemo::ContextArray& localValues ) const;

    virtual Scalar dotProduct( const Vector& other ) const;

    using Vector::prefetch; // prefetch() with no arguments

    virtual void prefetch( const hmemo::ContextPtr location ) const;

    virtual void wait() const;

    virtual void invert();

    virtual size_t getMemoryUsage() const;

    virtual void redistribute( DistributionPtr distribution );

    void writeToFile(
        const std::string& fileBaseName,
        const File::FileType fileType = File::BINARY,
        const File::DataType dataType = File::INTERNAL ) const;

protected:

    virtual void resizeImpl();

    using Vector::mContext;

    SCAI_LOG_DECL_STATIC_LOGGER( logger )

private    :

    void writeVectorToFormattedFile(const std::string& fileName) const;

    void writeVectorToBinaryFile(
                    const std::string& fileName,
                    const File::DataType outputType ) const;

    void writeVectorToXDRFile(
                    const std::string& fileName,
                    const File::DataType outputType ) const;

    void writeVectorDataToBinaryFile(
                    std::fstream& outFile,
                    const File::DataType outputType ) const;

    void readVectorHeader( const std::string& filename, File::FileType& fileType, long& dataTypeSize );

    void writeVectorHeader(
                    const std::string& fileName,
                    const File::FileType& fileType,
                    const long dataTypeSize ) const;

    void writeVectorToMMFile(
                    const std::string& filename,
                    const File::DataType& dataType ) const;

    void readVectorFromFormattedFile( const std::string& fileName );

    void readVectorFromBinaryFile(
                    const std::string& fileName,
                    const File::DataType dataType );

    void readVectorFromXDRFile(
                    const std::string& fileName,
                    const long dataTypeSizeHeader );

    void readVectorFromMMFile( const std::string& fileName );

    void readVectorDataFromBinaryFile(
                    std::fstream &inFile,
                    const File::DataType dataType );

    hmemo::LAMAArray<ValueType> mLocalValues; //!< my local values of vector

    mutable hmemo::LAMAArray<ValueType> mHaloValues;//!< my halo values of vector

public:

    // static methods, variables to register create routine in Vector factory of base class.

    static Vector* create();

    // key for factory 

    static std::pair<VectorKind, common::ScalarType> createValue();
};

/* ------------------------------------------------------------------------- */

template<typename ValueType>
template<typename OtherValueType>
DenseVector<ValueType>::DenseVector( const IndexType size, const OtherValueType* values, hmemo::ContextPtr context )
                : Vector( size, context )
{
    // use LAMA array reference to avoid copy of the raw data

    hmemo::LAMAArrayRef<OtherValueType> valuesArrayRef( size, values );

    // use mContext instead of context to avoid NULL pointer

    LAMAArrayUtils::assign( mLocalValues, valuesArrayRef, mContext );

    // Halo is not used yet

}

} /* end namespace lama */

} /* end namespace scai */
