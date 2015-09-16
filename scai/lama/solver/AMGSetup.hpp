/**
 * @file AMGSetup.hpp
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
 * @brief AMGSetup.hpp
 * @author Jiri Kraus
 * @date 28.10.2011
 * @since 1.0.0
 */

#pragma once

// for dll_import
#include <scai/common/config.hpp>

// local library
#include <scai/lama/matrix/Matrix.hpp>

#include <scai/lama/solver/Solver.hpp>

// internal scai libraries
#include <scai/common/Factory.hpp>

namespace scai
{

namespace lama
{

typedef common::shared_ptr<class AMGSetup> AMGSetupPtr;

/**
 * @brief The class AMGSetup should describe the Interace to an AMG Setup.
 *
 * @todo The current Interface of AMGSetup is just for evaluation so this should be changed to meet all requirements.
 *       (e.g. Pre and Post Smoothing)
 */
class COMMON_DLL_IMPORTEXPORT AMGSetup : 

    public scai::common::Factory<std::string, AMGSetup*>

{
public:

    AMGSetup();

    virtual ~AMGSetup();

    virtual void initialize( const Matrix& coefficients ) = 0;

    virtual Solver& getCoarseLevelSolver() = 0;

    virtual unsigned int getNumLevels() = 0;

    virtual Solver& getSmoother( const unsigned int level ) = 0;

    virtual const Matrix& getGalerkin( const unsigned int level ) = 0;

    virtual const Matrix& getRestriction( const unsigned int level ) = 0;

    virtual const Matrix& getInterpolation( const unsigned int level ) = 0;

    virtual Vector& getSolutionVector( const unsigned int level ) = 0;

    virtual Vector& getRhsVector( const unsigned int level ) = 0;

    virtual Vector& getTmpResVector( const unsigned int level ) = 0;

    virtual std::string getCouplingPredicateInfo() const = 0;

    virtual std::string getColoringInfo() const = 0;

    virtual std::string getInterpolationInfo() const = 0;

    virtual std::string getSmootherInfo() const = 0;

    virtual std::string getCoarseLevelSolverInfo() const = 0;

    virtual void setMaxLevels( const unsigned int level ) = 0;

    virtual void setMinVarsCoarseLevel( const unsigned int vars ) = 0;

    virtual void setHostOnlyLevel( IndexType hostOnlyLevel );

    virtual void setHostOnlyVars( IndexType hostOnlyVars );

    virtual void setReplicatedLevel( IndexType replicatedLevel );

    virtual void setCoarseLevelSolver( SolverPtr solver ) = 0;

    /**
     * @brief Sets smoother for all level
     */
    virtual void setSmoother( SolverPtr solver ) = 0;

protected:

    IndexType mHostOnlyLevel;

    IndexType mHostOnlyVars;

    IndexType mReplicatedLevel;
};

} /* end namespace lama */

} /* end namespace scai */
