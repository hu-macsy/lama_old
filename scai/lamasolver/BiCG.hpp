/**
 * @file BiCG.hpp
 *
 * @license
 * Copyright (c) 2013
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
 * @brief BiCG.hpp
 * @author Lauretta Schubert
 * @date 03.07.2013
 * @since 1.1.0
 */

#pragma once

// base classes
#include <scai/lamasolver/Solver.hpp>
#include <scai/lamasolver/CG.hpp>

// internal scai libraries
#include <scai/lama/matrix/SparseMatrix.hpp>

namespace scai
{

namespace lama
{

/**
 * @brief The class BiCG represents a IterativeSolver which uses the krylov subspace BiCG method
 *        to solve a system of linear equations iteratively.
 */
class COMMON_DLL_IMPORTEXPORT BiCG:

    public CG,
    public Solver::Register<BiCG>

{
public:

    /**
     * @brief Creates a BiCG solver with a given ID.
     *
     * @param id The ID for the solver.
     */
    BiCG( const std::string& id );

    /**
     * @brief Create a BiCG solver with a given ID and a given logger.
     *
     * @param id        The ID of the solver.
     * @param logger    The logger which shall be used by the solver
     */
    BiCG( const std::string& id, LoggerPtr logger );

    /**
     * @brief Copy constructor that copies the status independent solver information
     */
    BiCG( const BiCG& other );

    virtual ~BiCG();

    virtual void initialize( const Matrix& coefficients );

    /**
     * @brief Copies the status independent solver informations to create a new instance of the same
     * type
     *
     * @return shared pointer of the copied solver
     */
    virtual SolverPtr copy();

    struct BiCGRuntime: CGRuntime
    {
        BiCGRuntime();
        virtual ~BiCGRuntime();

        common::shared_ptr<Matrix> mTransposeA;
        common::shared_ptr<Vector> mP2;
        common::shared_ptr<Vector> mQ2;
        common::shared_ptr<Vector> mZ2;
        Scalar mPScalar2;
        mutable common::shared_ptr<Vector> mResidual2;
    };

    const Vector& getResidual2() const;

    /**
     * @brief Returns the complete configuration of the derived class
     */
    virtual BiCGRuntime& getRuntime();

    /**
     * @brief Returns the complete configuration of the derived class
     */
    virtual const BiCGRuntime& getConstRuntime() const;

    /**
     * @brief returns value used for registration of this solver
     */
    static std::string createValue();

    /**
     * @brief create a new BiCG solver with the corresponding name
     *
     * This method is used as create routine in the Solver factory.
     */
    static Solver* create( const std::string name );

protected:

    virtual void iterate();

    /**
     *  @brief own implementation of Printable::writeAt
     */
    virtual void writeAt( std::ostream& stream ) const;

    SCAI_LOG_DECL_STATIC_LOGGER( logger )

    BiCGRuntime    mBiCGRuntime;
};

} /* end namespace lama */

} /* end namespace scai */
