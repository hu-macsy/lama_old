/**
 * @file Criterion.cpp
 *
 * @license
 * Copyright (c) 2009-2016
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the Library of Accelerated Math Applications (LAMA).
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
 * @endlicense
 *
 * @brief Criterion.cpp
 * @author Kai Buschulte
 * @date 25.07.2011
 */

// hpp
#include <scai/solver/criteria/Criterion.hpp>

// local library
#include <scai/solver/IterativeSolver.hpp>

// internal scai libraries
#include <scai/common/macros/assert.hpp>

namespace scai
{

namespace solver
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

/* unused: operator= always used on CriterionPtr
Criterion& Criterion::operator=( const Criterion& other )
{
    mLeftChild = other.mLeftChild;
    mRightChild = other.mRightChild;
    mOperation = other.mOperation;

    return *this;
}*/

bool Criterion::isSatisfied( const IterativeSolver& solver )
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
    else if ( hasLeftChild() || hasRightChild() )
    {
        if( !mModifier )
        {
            stream << "!";
        }

        stream << "(";

        if ( hasLeftChild() )
        {
            stream << *getLeftChild();
        }
        else
        {
            stream << *getRightChild();
        }

        stream << ")";
    }
    else //leaf
    {
        stream << "leaf<" << mModifier << ">";
    }
}

} /* end namespace solver */

} /* end namespace scai */
