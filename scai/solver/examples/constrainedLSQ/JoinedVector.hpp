/**
 * @file JoinedVector.hpp
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
 * @brief Vector class that stands for the concatenation of two arbitrary vectors
 * @author Thomas Brandes
 * @date 27.07.2017
 */

#pragma once

#include <scai/lama/Vector.hpp>
#include <scai/lama/Scalar.hpp>

#include "JoinedDistribution.hpp"

namespace scai

{

namespace lama

{

/** An object of this class encapsulates a joined vector [ V1 ; V2 ] 
 */

class JoinedVector : public Vector
{

public:

    JoinedVector ( Vector& v1, Vector& v2 ) :

        Vector( dmemo::DistributionPtr( new dmemo::JoinedDistribution( v1.getDistributionPtr(), v2.getDistributionPtr() ) ) ),

        mV1( v1 ),
        mV2( v2 )

    {
        SCAI_ASSERT_EQ_ERROR( v1.getValueType(), v2.getValueType(), "joined vectors must have same type" )
    }

    JoinedVector ( Vector* v1, Vector* v2 ) :

        Vector( dmemo::DistributionPtr( new dmemo::JoinedDistribution( v1->getDistributionPtr(), v2->getDistributionPtr() ) ) ),

        mV1( *v1 ),
        mV2( *v2 )

    {
        SCAI_ASSERT_EQ_ERROR( v1->getValueType(), v2->getValueType(), "joined vectors must have same type" )
        pV1.reset( v1 );
        pV2.reset( v2 );
    }

    ~JoinedVector()
    {
    }

    Vector& first()
    {
        return mV1;
    }

    Vector& second()
    {
        return mV2;
    }

    const Vector& first() const
    {
        return mV1;
    }

    const Vector& second() const
    {
        return mV2;
    }

 	virtual bool isConsistent() const
    {
        COMMON_THROWEXCEPTION( "unsupported" )
    }

	virtual _Vector::VectorKind getVectorKind() const
    {
        return Vector::JOINED;
    }

	virtual void setDenseValues(const scai::hmemo::_HArray&)
    {
        COMMON_THROWEXCEPTION( "unsupported" )
    }

	virtual void setSparseValues(const scai::hmemo::HArray<int>&, const scai::hmemo::_HArray&, Scalar)
    {
        COMMON_THROWEXCEPTION( "unsupported" )
    }

	virtual void clearValues()
    {
        COMMON_THROWEXCEPTION( "unsupported" )
    }

	virtual void setRandom(scai::dmemo::DistributionPtr, float)
    {
        COMMON_THROWEXCEPTION( "unsupported" )
    }

	virtual scai::common::scalar::ScalarType getValueType() const
    {
        COMMON_THROWEXCEPTION( "unsupported" )
    }

	virtual Scalar getValue(IndexType) const
    {
        COMMON_THROWEXCEPTION( "unsupported" )
    }

	virtual void setValue(IndexType, Scalar)
    {
        COMMON_THROWEXCEPTION( "unsupported" )
    }

	virtual Scalar min() const
    {
        COMMON_THROWEXCEPTION( "unsupported" )
    }

	virtual Scalar max() const
    {
        COMMON_THROWEXCEPTION( "unsupported" )
    }

	virtual Scalar sum() const
    {
        COMMON_THROWEXCEPTION( "unsupported" )
    }

	virtual Scalar l1Norm() const
    {
        COMMON_THROWEXCEPTION( "unsupported" )
    }

	virtual Scalar l2Norm() const
    {
        Scalar res = dotProduct( *this );
        return lama::sqrt( res );
    }

	virtual Scalar maxNorm() const
    {
        COMMON_THROWEXCEPTION( "unsupported" )
    }

	virtual Scalar maxDiffNorm(const Vector&) const
    {
        COMMON_THROWEXCEPTION( "unsupported" )
    }

	virtual VectorCreateKeyType getCreateValue() const
    {
        COMMON_THROWEXCEPTION( "unsupported" )
    }

	virtual Vector* newVector() const
    {
        COMMON_THROWEXCEPTION( "unsupported" )
    }

	virtual Vector* copy() const
    {
        COMMON_THROWEXCEPTION( "unsupported" )
    }

	virtual void swap(Vector&)
    {
        COMMON_THROWEXCEPTION( "unsupported" )
    }

	virtual void buildLocalValues(scai::hmemo::_HArray&, scai::common::BinaryOp, scai::hmemo::ContextPtr) const
    {
        COMMON_THROWEXCEPTION( "unsupported" )
    }

	virtual void gatherLocalValues(scai::hmemo::_HArray&, const scai::hmemo::HArray<int>&, scai::common::BinaryOp, scai::hmemo::ContextPtr) const
    {
        COMMON_THROWEXCEPTION( "unsupported" )
    }

	virtual void assign( Scalar s )
    {
        mV1.assign( s );
        mV2.assign( s );
    }

	virtual void vectorPlusVector( const Scalar& alpha, const Vector& v1, const Scalar& beta, const Vector& v2 )
    {
        SCAI_ASSERT_EQ_ERROR( v1.getDistribution(), getDistribution(), "distribution mismatch" );
        SCAI_ASSERT_EQ_ERROR( v2.getDistribution(), getDistribution(), "distribution mismatch" );

        const JoinedVector& jV1 = reinterpret_cast<const JoinedVector&>( v1 );
        const JoinedVector& jV2 = reinterpret_cast<const JoinedVector&>( v2 );
    
        mV1.vectorPlusVector( alpha, jV1.first(), beta, jV2.first() );
        mV2.vectorPlusVector( alpha, jV1.second(), beta, jV2.second() );
    }

	virtual void vectorTimesVector(const Scalar&, const Vector&, const Vector&)
    {
        COMMON_THROWEXCEPTION( "unsupported" )
    }

	virtual void vectorPlusScalar(const Scalar&, const Vector&, const Scalar&)
    {
        COMMON_THROWEXCEPTION( "unsupported" )
    }

	virtual Scalar dotProduct( const Vector& other ) const
    {
        SCAI_ASSERT_EQ_ERROR( other.getDistribution(), getDistribution(), "distribution mismatch" );
        
        const JoinedVector& jOther = reinterpret_cast<const JoinedVector&>( other );

        Scalar s1 = mV1.dotProduct( jOther.first() );
        Scalar s2 = mV2.dotProduct( jOther.second() );

        return mV1.dotProduct( jOther.first() ) + mV2.dotProduct( jOther.second() );
    }

	virtual void setVector(const Vector&, scai::common::BinaryOp, bool)
    {
        COMMON_THROWEXCEPTION( "unsupported" )
    }

	virtual void setScalar(Scalar, scai::common::BinaryOp, bool)
    {
        COMMON_THROWEXCEPTION( "unsupported" )
    }

	virtual void applyUnary(scai::common::unary::UnaryOp)
    {
        COMMON_THROWEXCEPTION( "unsupported" )
    }

	virtual bool all(scai::common::CompareOp, Scalar) const
    {
        COMMON_THROWEXCEPTION( "unsupported" )
    }

	virtual bool all(scai::common::CompareOp, const Vector&) const
    {
        COMMON_THROWEXCEPTION( "unsupported" )
    }

	virtual void prefetch(scai::hmemo::ContextPtr) const
    {
        COMMON_THROWEXCEPTION( "unsupported" )
    }

	virtual void wait() const
    {
        COMMON_THROWEXCEPTION( "unsupported" )
    }

	virtual size_t getMemoryUsage() const
    {
        COMMON_THROWEXCEPTION( "unsupported" )
    }

	virtual void allocate(scai::dmemo::DistributionPtr)
    {
        COMMON_THROWEXCEPTION( "unsupported" )
    }

	virtual void allocate(IndexType)
    {
        COMMON_THROWEXCEPTION( "unsupported" )
    }

	virtual void redistribute(scai::dmemo::DistributionPtr)
    {
        COMMON_THROWEXCEPTION( "unsupported" )
    }

	virtual void writeLocalToFile(const std::string&, const std::string&, scai::common::scalar::ScalarType, FileIO::FileMode) const
    {
        COMMON_THROWEXCEPTION( "unsupported" )
    }

	virtual IndexType readLocalFromFile(const std::string&, IndexType, IndexType)
    {
        COMMON_THROWEXCEPTION( "unsupported" )
    }

    /** Override default implementation of Printable::writeAt */

    virtual void writeAt( std::ostream& stream ) const
    {
        stream << "JoinedVector( " << mV1 << ", " << mV2 << " )";
    }

private:

    Vector& mV1;
    Vector& mV2;

    // joined vector might be owner of the individual parts

    _VectorPtr pV1;
    _VectorPtr pV2;
};

}

}
