/**
 * @file ResidualStagnation.hpp
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
 * @brief ResidualStagnation.hpp
 * @author Kai Buschulte
 * @date 25.07.2011
 * $Id$
 */

#ifndef LAMA_RESIDUALSTAGNATION_HPP_
#define LAMA_RESIDUALSTAGNATION_HPP_

// for dll_import
#include <lama/config.hpp>

// base classes
#include <lama/solver/criteria/Criterion.hpp>

// others
#include <lama/LAMATypes.hpp>

#include <lama/norm/Norm.hpp>

namespace lama
{

class IterativeSolver;

/**
 * @brief Stagnation stopping criterion.
 *
 * @param IterativeSolver The type solver in which this criterion is used.
 *        Most commonly this is the abstract sblas::IterativeSolver type.
 */
class LAMA_DLL_IMPORTEXPORT ResidualStagnation: public Criterion
{

public:

    /**
     * @brief Creates a new Stagnation based stopping criterion.
     * TODO[doxy] param in or out?
     *
     * @param      norm The norm which shall be used to check the criterion.
     *             Has to be compatible with the Matrix/Vector types of
     *             the solver.
     */
    ResidualStagnation( NormPtr norm );

    /**
     * @brief Creates a new Stagnation based stopping criterion.
     * TODO[doxy] params in or out?
     *
     * @param      norm The norm which shall be used to check the criterion.
     *             Has to be compatible with the Matrix/Vector types of
     *             the solver.
     * @param      lookback The amount of residual-norm-calculation results used
     *             for the criterion check.
     * @param      precision TODO[doxy] Complete Description.
     */
    ResidualStagnation( NormPtr norm, IndexType lookback, Scalar precision );

    /**
     * @brief Creates a copy of the passed ResidualStagnation object.
     *
     * @param[in] other   ResidualStagnation object to be copied.
     */
    ResidualStagnation( const ResidualStagnation &other );

    /** Destructor. */

    virtual ~ResidualStagnation();

    /**
     * @brief TODO[doxy] Complete Description.
     *
     * @return   TODO[doxy] Complete Description.
     */
    virtual Criterion* copy() const;

    /**
     * @brief Checks if the criterion is satisfied.
     *
     * The Stagnation criterion needs a certain amount (specified by
     * m_lookback) to be able to perform the check. If this method is
     * called and not enough data has been gathered, false is returned.
     * Anyhow the norm of the residual is calculated and stored. This
     * process continues, until enough data has been gathered.
     *
     * @param[in] solver   The solver which calls this method. Needed for
     *                     residual computation.
     * @return   Boolean value which determines if the criterion is
     *           satisfied or not.
     */
    virtual bool isSatisfied( const IterativeSolver& solver );

    /**
     * @brief Gets the amount of calculations used for the criterion check.
     *
     * @return   Amount of calculations.
     */
    IndexType getLookback() const;

    /**
     * @brief Gets the Norm which the criterion uses for the check.
     *
     * @return   The Norm used.
     */
    const NormPtr getNorm() const;

    /**
     * @brief Returns the precision for floating point number calculations related to the Norm.
     *
     * @return   The precision given in a Scalar object.
     */
    const Scalar getPrecision() const;

    /**
     * @brief Sets the amount of calculations used for the criterion check.
     *
     * @param[in] lookback   TODO[doxy] Complete Description.
     */
    void setLookback( const IndexType lookback );

    /**
     * @brief Sets the precision for floating point number calculations related to the Norm.
     *
     * @param[in] precision   The precision given in a Scalar object.
     */
    void setPrecision( const Scalar precision );

    virtual void writeAt( std::ostream& stream ) const;

private:

    /**
     * @brief The norm used in the residual-norm-calculation.
     */
    const NormPtr mNorm;

    /**
     * @brief The amount of calculations used for the criterion check.
     */
    IndexType mLookback;

    /**
     * @brief Stores the results of the residual-norm-calculations.
     */
    std::vector<Scalar> mLastResidualNorms;

    /**
     * @brief Index. Needed for circular use of the vector.
     */
    IndexType mNextEntry;

    /**
     * @brief Used to check, if enough data have been gathered to check the
     *        criterion.
     */
    bool mEntriesReady;

    /**
     * @brief The precision used for the stagnation check.
     */
    Scalar mPrecision;
};

} // namespace lama

#endif // LAMA_RESIDUALSTAGNATION_HPP_
