/**
 * @file SpecializedJacobi.hpp
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
 * @brief SpecializedJacobi.hpp
 * @author Matthias Makulla
 * @date 06.04.2011
 * $Id$
 */

#ifndef LAMA_SPECIALIZEDJACOBI_HPP_
#define LAMA_SPECIALIZEDJACOBI_HPP_

// for dll_import
#include <lama/config.hpp>

// base classes
#include <lama/solver/OmegaSolver.hpp>

// others
#include <lama/matrix/SparseMatrix.hpp>

namespace lama
{

class LAMA_DLL_IMPORTEXPORT SpecializedJacobi: public OmegaSolver
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

        //TODO: LAMAArray?
        std::auto_ptr<Vector> mOldSolution;
        SolutionProxy mProxyOldSolution;
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

protected:
    SpecializedJacobiRuntime mSpecializedJacobiRuntime;
private:
    template<typename ValueType>
    void iterateTyped( const SparseMatrix<ValueType>& );
};

} // namespace lama

#endif // LAMA_SPECIALIZEDJACOBI_HPP_
