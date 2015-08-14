/**
 * @file Criterion.hpp
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
 * @brief Criterion.hpp
 * @author Kai Buschulte
 * @date 21.07.2011
 * @since 1.0.0
 */

#ifndef LAMA_BOOLEANCONDITION_HPP_
#define LAMA_BOOLEANCONDITION_HPP_

// for dll_import
#include <scai/common/config.hpp>

// base classes
#include <scai/common/Printable.hpp>

// others
#include <scai/lama/LAMATypes.hpp>

// logging
#include <scai/logging.hpp>

// boost
#include <scai/common/shared_ptr.hpp>

#include <ostream>

namespace lama
{

class IterativeSolver;

class Criterion;

typedef common::shared_ptr<Criterion> CriterionPtr;

/**
 * @brief The class Criterion is the base class for all stopping criterions
 *        that we can use on a Solver.
 *
 * Criterion can be composed with the logical operators &&, ||, and !.
 */
class COMMON_DLL_IMPORTEXPORT Criterion: public Printable
{
public:

    /**
     * @brief Defines the operators which can be used to connect Criterions with.
     */
    enum BooleanOperator
    {
        AND, //!< stands for logical and composition
        OR //!< stands for logical or composition
    };

    /**
     * @brief Creates a Criterion which evaluates to true.
     */
    Criterion();

    /**
     * @brief Creates a Criterion which evaluates to the passed boolean.
     *
     * @param[in] boolean   the boolean to set the result of this boolean condition to.
     *
     * Note: implicit type conversions of bool argument are disabled
     */
    explicit Criterion( const bool boolean );

    /**
     * @brief Creates a copy of the passed boolean condition.
     *
     * @param[in] other the boolean condition to copy.
     */
    Criterion( const Criterion& other );

    /**
     * @brief Creates a copy of the passed Criterion, which is altered by the
     *        passed modifier.
     *
     * Creates a copy of the passed Criterion, which is altered by the
     * passed modifier. The new Criterion behaves like
     * \code
     *     other.isSatisfied() != modifier
     * \endcode
     * which is equivalent to  other.isSatisfied() XOR modifier.
     *
     * @param[in] other     the boolean condition to take a altered copy from
     * @param[in] modifier  the bool value to modify the new Criterion with.
     */
    Criterion( const CriterionPtr other, const bool modifier );

    /**
     * @brief Creates a new Criterion which the logical composition of
     *        the passed Criterions.
     *
     * @param[in] leftChild     the left operand of the composition
     * @param[in] rightChild    the right operand of the composition
     * @param[in] operation     the logical operation to connect leftChild and
     *                          rightChild with.
     */
    Criterion( const CriterionPtr leftChild, const CriterionPtr rightChild, BooleanOperator operation );

    /** Destructor. */

    virtual ~Criterion();

    /**
     * @brief Checks if this Criterion is true for the state of the passed
     *        solver.
     *
     * @param[in] solver    the solver to check the state of.
     * @return              if this Criterion is true for the state of solver.
     */
    virtual bool isSatisfied( const IterativeSolver& solver );

    const CriterionPtr getLeftChild() const;

    const CriterionPtr getRightChild() const;

    BooleanOperator getOperation() const;

    void setLeftChild( const CriterionPtr leftChild );

    void setRightChild( const CriterionPtr rightChild );

    void setOperation( const BooleanOperator operation );

    /**
     * @brief Returns true if the boolean condition has a left child
     * Do not use getLeftChild() to compare pointers. Logical operations are overwritten
     *
     * return Returns true if the boolean condition has a left child. If pointer is not set: false
     */
    bool hasLeftChild() const;

    /**
     * @brief Returns true if the boolean condition has a right child
     * Do not use getRightChild() to compare pointers. Logical operations are overwritten
     *
     * return Returns true if the boolean condition has a right child. If pointer is not set: false
     */
    bool hasRightChild() const;

    Criterion& operator =( const Criterion& other );

    virtual void writeAt( std::ostream& stream ) const;

protected:

    SCAI_LOG_DECL_STATIC_LOGGER( logger )

private    :
    CriterionPtr mLeftChild;
    CriterionPtr mRightChild;
    BooleanOperator mOperation;
    const bool mModifier; // true: negate, false: do nothing
};

inline CriterionPtr operator||( CriterionPtr a, CriterionPtr b )
{
    return CriterionPtr( new Criterion( a, b, Criterion::OR ) );
}

inline CriterionPtr operator&&( CriterionPtr a, CriterionPtr b )
{
    return CriterionPtr( new Criterion( a, b, Criterion::AND ) );
}

} // namespace lama

#endif // LAMA_BOOLEANCONDITION_HPP_
