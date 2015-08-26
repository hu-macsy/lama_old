/**
 * @file SingleGridSetup.hpp
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
 * @brief SingleGridSetup.hpp
 * @author Jiri Kraus
 * @date 27.10.2011
 * @since 1.0.0
 */

#pragma once

// for dll_import
#include <scai/common/config.hpp>

// base classes
#include <scai/lama/solver/AMGSetup.hpp>

// others
#include <scai/lama/solver/Solver.hpp>

namespace scai
{

namespace lama
{

class SingleGridSetup: 

    public AMGSetup,
    public AMGSetup::Register<SingleGridSetup>    // register at factory

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
    // just a dummy function
    virtual void setMinVarsCoarseLevel( const unsigned int )
    {
    }

    // Get the key used for registration in factory

    static std::string createValue();

    // Create routine used for factory 

    static AMGSetup* create();

private:

    SCAI_LOG_DECL_STATIC_LOGGER( logger )

    SolverPtr    mSolver;
    common::unique_ptr<Matrix> mIdentity;

    common::unique_ptr<Vector> mSolutionVector;
    common::unique_ptr<Vector> mRhsVector;
    common::unique_ptr<Vector> mTmpResVector;

};

} /* end namespace lama */

} /* end namespace scai */
