/**
 * @file IterationCount.hpp
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
 * @brief IterationCount.hpp
 * @author Kai Buschulte
 * @date 21.07.2011
 * @since 1.0.0
 */

#pragma once

// for dll_import
#include <scai/common/config.hpp>

// base classes
#include <scai/solver/criteria/Criterion.hpp>

namespace scai
{

namespace solver
{

class IterativeSolver;

/**
 * @brief IterationCount is a stopping criterion of a solver which checks the
 *        number of iterations a solver has executed.
 *
 * IterationCount is a stopping criterion of a solver which checks the number of
 * iterations a solver has executed. IterationCount is either true if the number
 * of iterations is larger or smaller than a configured number of iterations.
 * Depending on the configured IterationCheckMode.
 */
class COMMON_DLL_IMPORTEXPORT IterationCount: public Criterion
{
public:
    /**
     * @brief Creates a IterationCount with iterationExtream 1.
     */
    IterationCount();

    /**
     * @brief Creates a IterationCount with the passed iterationExtrema.
     *
     * @param[in] iterationExtrema   the number of iterations a solver should execute at most or minimal
     */
    IterationCount( const IndexType iterationExtrema );

    /**
     * @brief Creates a copy of the passed IterationCount object.
     *
     * @param[in] other   IterationCount object to be copied.
     */
    IterationCount( const IterationCount &other );

    /** Destructor. */

    virtual ~IterationCount();

    /**
     * @brief TODO[doxy] Complete Description.
     *
     * @param[in] solver   TODO[doxy] Complete Description.
     * @return             TODO[doxy] Complete Description.
     */
    virtual bool isSatisfied( const IterativeSolver& solver );

    /**
     * @brief Getter of the iteration extrema.
     *
     * @return   the iteration extrema.
     */
    IndexType getIterationExtrema() const;

    /**
     * @brief Setter of the iteration extrema.
     *
     * @param[in] iterationExtrema   the new iteration extrema.
     */
    void setIterationExtrema( IndexType iterationExtrema );

    virtual void writeAt( std::ostream& stream ) const;

protected:

    SCAI_LOG_DECL_STATIC_LOGGER( logger )

private    :

    IndexType mIterationExtrema;
};

} /* end namespace solver */

} /* end namespace scai */
