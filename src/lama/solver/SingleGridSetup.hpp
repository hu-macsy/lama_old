/**
 * @file SingleGridSetup.hpp
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
 * @brief SingleGridSetup.hpp
 * @author Jiri Kraus
 * @date 27.10.2011
 * $Id$
 */
#ifndef LAMA_SINGLEGRIDSETUP_HPP_
#define LAMA_SINGLEGRIDSETUP_HPP_

// for dll_import
#include <lama/config.hpp>

// base classes
#include <lama/solver/AMGSetup.hpp>

// others
#include <lama/solver/Solver.hpp>

namespace lama
{

class SingleGridSetup: public lama::AMGSetup
{
public:
    SingleGridSetup();
    virtual ~SingleGridSetup();

    virtual void initialize( const Matrix& coefficients );

    virtual Solver& getCoarseLevelSolver();

    virtual unsigned int getNumLevels();

    virtual Solver& getSmoother( const unsigned int level );

    virtual const Matrix& getGalerkin( const unsigned int level );

    virtual const Matrix& getRestriction( const unsigned int level );

    virtual const Matrix& getInterpolation( const unsigned int level );

    virtual Vector& getSolutionVector( const unsigned int level );

    virtual Vector& getRhsVector( const unsigned int level );

    virtual Vector& getTmpResVector( const unsigned int level );

    virtual std::string getCouplingPredicateInfo() const;

    virtual std::string getColoringInfo() const;

    virtual std::string getInterpolationInfo() const;

    virtual std::string getSmootherInfo() const;

    virtual std::string getCoarseLevelSolverInfo() const;

    virtual void setCoarseLevelSolver( SolverPtr solver );

    /**
     * @brief Usually sets smoother for all level, for this case overwrites the coarse level solver
     */
    virtual void setSmoother( SolverPtr solver );

    // just a dummy function
    virtual void setMaxLevels( const unsigned int )
    {
    }
    ;
    // just a dummy function
    virtual void setMinVarsCoarseLevel( const unsigned int )
    {
    }
    ;

private:

    LAMA_LOG_DECL_STATIC_LOGGER( logger )

    SolverPtr mSolver;
    std::auto_ptr<Matrix> mIdentity;

    std::auto_ptr<Vector> mSolutionVector;
    std::auto_ptr<Vector> mRhsVector;
    std::auto_ptr<Vector> mTmpResVector;

};

}

#endif /* LAMA_SINGLEGRIDSETUP_HPP_ */
