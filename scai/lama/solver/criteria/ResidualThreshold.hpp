/**
 * @file ResidualThreshold.hpp
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
 * @brief ResidualThreshold.hpp
 * @author Kai Buschulte
 * @date 21.07.2011
 * @since 1.0.0
 */

#ifndef LAMA_RESIDUALTHRESHOLD_HPP_
#define LAMA_RESIDUALTHRESHOLD_HPP_

// for dll_import
#include <scai/common/config.hpp>

// base classes
#include <scai/lama/solver/criteria/Criterion.hpp>

// others
#include <scai/lama/Scalar.hpp>

#include <scai/lama/norm/Norm.hpp>

namespace lama
{

/**
 * @brief ResidualThreshold is a stopping criterion of a solver which checks the
 *        residue of a solver for certain criteria.
 *
 * ResidualThreshold is a stopping criterion of a solver which checks the residue
 * of a solver to fall below a configured threshold. The threshold is either checked
 * absolute or relative to the initial residue.
 *
 */
class COMMON_DLL_IMPORTEXPORT ResidualThreshold: public Criterion
{
public:
    /**
     * @brief The possible check modes for the residue threshold
     */
    enum ResidualThresholdCheckMode
    {
        Absolute, //!< stands for an absolute residue reduction
        Relative //!< stands for an relative residue reduction
    };

    /*
     * @brief Empty constructor first. Configured by setter methods.
     */
    ResidualThreshold();

    /**
     * @brief Creates a RedisualThreshold stopping criteria from the passed norm.
     *
     * @param[in] norm      the norm to use for residue calculation
     */
    ResidualThreshold( const NormPtr norm );

    /**
     * @brief Creates a RedisualThreshold stopping criteria from the passed parameters.
     *
     * @param[in] norm      the norm to use for residue calculation
     * @param[in] precision the desired precision
     * @param[in] checkMode if the residue should be check for an absolute or a relative
     *                      reduction.
     */
    ResidualThreshold( const NormPtr norm, Scalar precision, ResidualThresholdCheckMode checkMode );

    /**
     * @brief Creates a copy of the passed ResidualThreshold.
     *
     * @param[in] other the ResidualThreshold to take a copy from.
     */
    ResidualThreshold( const ResidualThreshold& other );

    /**
     * @brief Destroys this ResidualThreshold.
     */
    virtual ~ResidualThreshold();

    virtual Criterion* copy() const;

    bool isSatisfied( const IterativeSolver& solver );

    ResidualThresholdCheckMode getCheckMode() const;
    IndexType getIterationExtrema() const;
    Scalar getFirstNormResult() const;
    const NormPtr getNorm() const;
    Scalar getPrecision() const;

    void setCheckMode( ResidualThresholdCheckMode mCheckMode );
    void setFirstNormResult( Scalar firstNormResult );
    void setPrecision( Scalar precision );

    virtual void writeAt( std::ostream& stream ) const;

protected:
    LAMA_LOG_USING( Criterion::logger );
private:
    const NormPtr mNorm;
    ResidualThresholdCheckMode mCheckMode;
    Scalar mPrecision;
    Scalar mFirstNormResult;
};

}
//namespace lama

#endif // LAMA_RESIDUALTHRESHOLD_HPP_
