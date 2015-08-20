/**
 * @file Criterion.cpp
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
 * @brief Criterion.cpp
 * @author Kai Buschulte
 * @date 25.07.2011
 * @since 1.0.0
 */

// hpp
#include <scai/lama/solver/criteria/Criterion.hpp>

// others
#include <scai/lama/solver/IterativeSolver.hpp>

// assert
#include <scai/common/SCAIAssert.hpp>

namespace scai
{

namespace lama
{

SCAI_LOG_DEF_LOGGER( Criterion::logger, "Criterion" );

Criterion::Criterion()
    : mOperation( AND ), mModifier( true )
{
}

Criterion::Criterion( const bool boolean )
    : mOperation( AND ), mModifier( boolean )
{
}

Criterion::Criterion( const Criterion& other )
    : Printable( other ), mLeftChild( other.mLeftChild ), mRightChild( other.mRightChild ), mOperation(
          other.mOperation ), mModifier( false )
{
}

Criterion::Criterion( const CriterionPtr other, const bool modifier )
    : mLeftChild( other->mLeftChild ), mRightChild( other->mRightChild ), mOperation( other->mOperation ), mModifier(
          modifier == other->mModifier )
{
}

Criterion::Criterion( const CriterionPtr leftChild, const CriterionPtr rightChild, BooleanOperator operation )
    : mLeftChild( leftChild ), mRightChild( rightChild ), mOperation( operation ), mModifier( false )
{
    SCAI_ASSERT_DEBUG( leftChild.get() != NULL, "Left child is a NULL pointer." );
    SCAI_ASSERT_DEBUG( rightChild.get() != NULL, "Right child is a NULL pointer." );
}

Criterion::~Criterion()
{
}

Criterion& Criterion::operator=( const Criterion& other )
{
    mLeftChild = other.mLeftChild;
    mRightChild = other.mRightChild;
    mOperation = other.mOperation;

    return *this;
}

bool Criterion::isSatisfied( const scai::lama::IterativeSolver& solver )
{
    SCAI_LOG_INFO( logger, "isSatisfied: " << *this );

    bool satisfied = true;

    if( !hasLeftChild() && !hasRightChild() )
    {
        return mModifier;
    }

    if( !hasLeftChild() )
    {
        satisfied = mRightChild->isSatisfied( solver );
    }
    else if( !hasRightChild() )
    {
        satisfied = mLeftChild->isSatisfied( solver );
    }
    else if( mOperation == AND )
    {
        satisfied = mLeftChild->isSatisfied( solver ) && mRightChild->isSatisfied( solver );
    }
    else
    {
        satisfied = mLeftChild->isSatisfied( solver ) || mRightChild->isSatisfied( solver );
    }

    return mModifier != satisfied; //XOR
}

const CriterionPtr Criterion::getLeftChild() const
{
    return mLeftChild;
}

const CriterionPtr Criterion::getRightChild() const
{
    return mRightChild;
}

Criterion::BooleanOperator Criterion::getOperation() const
{
    return mOperation;
}

void Criterion::setLeftChild( const CriterionPtr leftChild )
{
    mLeftChild = leftChild;
}

void Criterion::setRightChild( const CriterionPtr rightChild )
{
    mRightChild = rightChild;
}

bool Criterion::hasLeftChild() const
{
    return mLeftChild.get();
}

bool Criterion::hasRightChild() const
{
    return mRightChild.get();
}

void Criterion::setOperation( const Criterion::BooleanOperator operation )
{
    mOperation = operation;
}

void Criterion::writeAt( std::ostream& stream ) const
{
    if( hasLeftChild() && hasRightChild() ) //boolean operation
    {
        if( !mModifier )
        {
            stream << "!";
        }

        stream << "(" << *getLeftChild();

        if( getOperation() == Criterion::AND )
        {
            stream << " && ";
        }
        else
        {
            stream << " || ";
        }

        stream << *getRightChild() << ")";
    }
    else //leaf
    {
        stream << "leaf<" << mModifier << ">";
    }
}

} /* end namespace lama */

} /* end namespace scai */
