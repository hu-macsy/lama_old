/**
 * @file SpecializedJacobi.hpp
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
 * @brief SpecializedJacobi.hpp
 * @author Matthias Makulla
 * @date 06.04.2011
 * @since 1.0.0
 */

#pragma once

// for dll_import
#include <scai/common/config.hpp>

// base classes
#include <scai/lama/solver/Solver.hpp>
#include <scai/lama/solver/OmegaSolver.hpp>

// local library
#include <scai/lama/matrix/SparseMatrix.hpp>

namespace scai
{

namespace lama
{

class COMMON_DLL_IMPORTEXPORT SpecializedJacobi:
	public OmegaSolver,
	public Solver::Register<SpecializedJacobi>
{
public:
    SpecializedJacobi( const std::string& id );
    SpecializedJacobi( const std::string& id, Scalar omega );
    SpecializedJacobi( const std::string& id, LoggerPtr logger );
    SpecializedJacobi( const std::string& id, Scalar omega, LoggerPtr logger );
    SpecializedJacobi( const SpecializedJacobi& other );

    virtual ~SpecializedJacobi();

    virtual void initialize( const Matrix& coefficients );

    virtual void solve( Vector& rhs, const Vector& solution );

    virtual void solveInit( Vector& solution, const Vector& rhs );

    virtual void solveFinalize();

    void iterate();

    struct SpecializedJacobiRuntime: OmegaSolverRuntime
    {
        SpecializedJacobiRuntime();
        virtual ~SpecializedJacobiRuntime();

        //TODO: HArray?
        common::shared_ptr<Vector> mOldSolution;
        SolutionProxy mProxyOldSolution;
        common::shared_ptr<hmemo::ContextArray> mDiagonal;
    };

    /**
     * @brief Returns the complete configuration of the derived class
     */
    virtual SpecializedJacobiRuntime& getRuntime();

    /**
     * @brief Returns the complete const configuration of the derived class
     */
    virtual const SpecializedJacobiRuntime& getConstRuntime() const;

    /**
     * @brief Copies the status independent solver informations to create a new instance of the same
     * type
     *
     * @return shared pointer of the copied solver
     */
    virtual SolverPtr copy();

    static std::string createValue();
    static Solver* create( const std::string name );

protected:
    SpecializedJacobiRuntime mSpecializedJacobiRuntime;
private:
    template<typename ValueType>
    void iterateTyped( const SparseMatrix<ValueType>& );

    template<typename ValueType>
    void iterateSync(
        hmemo::HArray<ValueType>& solution,
        const SparseMatrix<ValueType>& coefficients,
        const hmemo::HArray<ValueType>& localOldSolution,
        hmemo::HArray<ValueType>& haloOldSolution,
        const hmemo::HArray<ValueType>& rhs,
        const ValueType omega );

    template<typename ValueType>
    void iterateAsync(
        hmemo::HArray<ValueType>& solution,
        const SparseMatrix<ValueType>& coefficients,
        const hmemo::HArray<ValueType>& localOldSolution,
        hmemo::HArray<ValueType>& haloOldSolution,
        const hmemo::HArray<ValueType>& rhs,
        const ValueType omega );
};

} /* end namespace lama */

} /* end namespace scai */
