/**
 * @file GMRES.hpp
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
 * @brief GMRES.hpp
 * @author Malte Förster
 * @date 10.04.2012
 * @since 1.0.0
 */
#ifndef LAMA_GMRES_HPP_
#define LAMA_GMRES_HPP_

// for dll_import
#include <scai/common/config.hpp>

// base classes
#include <scai/lama/solver/IterativeSolver.hpp>

// boost
#include <scai/common/unique_ptr.hpp>

namespace scai
{

namespace lama
{

/**
 * @brief The class GMRES represents a IterativeSolver which uses the krylov subspace GMRES method
 *        to solve a system of linear equations iteratively.
 */
class COMMON_DLL_IMPORTEXPORT GMRES: public IterativeSolver
{
public:

    /**
     * @brief Creates a GMRES solver with a given ID.
     *
     * @param id The ID for the solver.
     */
    GMRES( const std::string& id );

    /**
     * @brief Create a GMRES solver with a given ID and a given logger.
     *
     * @param id        The ID of the solver.
     * @param logger    The logger which shall be used by the solver
     */
    GMRES( const std::string& id, LoggerPtr logger );

    /**
     * @brief Copy constructor that copies the status independent solver information
     */
    GMRES( const GMRES& other );

    virtual ~GMRES();

    virtual void initialize( const Matrix& coefficients );

    void setKrylovDim( unsigned int krylovDim );

    /**
     * @brief Copies the status independent solver informations to create a new instance of the same
     * type
     *
     * @return shared pointer of the copied solver
     */
    virtual SolverPtr copy();

    struct GMRESRuntime: IterativeSolverRuntime
    {
        GMRESRuntime();
        virtual ~GMRESRuntime();

        // arrays to store rotations
        common::scoped_array<double> mCC;
        common::scoped_array<double> mSS;

        // array for Hessenberg equation
        // H*y=g
        common::scoped_array<double> mG;
        common::scoped_array<double> mY;

        // Hessenberg matrix
        // mH:  Upper triangular (columnwise)
        // mHd: diagonal band h(i+1,i)
        common::scoped_array<double> mH;
        common::scoped_array<double> mHd;

        // krylov space
        std::vector<Vector*> *mV;

        // temp-arrays
        Vector *mW;
        Vector *mT;

        // remember starting solution
        // only needed if x is modified within krylov loop
        Vector *mX0;
    };

    /**
     * @brief Returns the complete configuration of the derived class
     */
    virtual GMRESRuntime& getRuntime();

    /**
     * @brief Returns the complete configuration of the derived class
     */
    virtual const GMRESRuntime& getConstRuntime() const;

    double getAverageIterationTime() const;
    double getAveragePreconditionerTime() const;

protected:

    virtual void iterate();

    GMRESRuntime mGMRESRuntime;

    SCAI_LOG_DECL_STATIC_LOGGER( logger )

private    :

    void updateX( unsigned int i );

    // krylov dimension
    unsigned int mKrylovDim;

    double totalIterationTime;
    double totalPreconditionerTime;
};

} /* end namespace lama */

} /* end namespace scai */

#endif // LAMA_GMRES_HPP_
