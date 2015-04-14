/**
 * @file TrivialPreconditioner.hpp
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
 * @brief TrivialPreconditioner.hpp
 * @author Matthias Makulla
 * @date 06.04.2011
 * @since 1.0.0
 */

#ifndef LAMA_TRIVIALPRECONDITIONER_HPP_
#define LAMA_TRIVIALPRECONDITIONER_HPP_

// for dll_import
#include <lama/config.hpp>

// base classes
#include <lama/solver/Solver.hpp>

namespace lama
{

class LAMA_DLL_IMPORTEXPORT TrivialPreconditioner: public Solver
{
public:
    TrivialPreconditioner( const std::string& id );
    TrivialPreconditioner( const std::string& id, LoggerPtr logger );
    TrivialPreconditioner( const TrivialPreconditioner& other );

    virtual ~TrivialPreconditioner();

    void initialize( const Matrix& coefficients );

    virtual void solveImpl();

    virtual SolverPtr copy();

    struct TrivialPreconditionerRuntime: SolverRuntime
    {
        TrivialPreconditionerRuntime();
        virtual ~TrivialPreconditionerRuntime();
    };

    virtual TrivialPreconditionerRuntime& getRuntime();

    virtual const TrivialPreconditionerRuntime& getConstRuntime() const;

protected:
    TrivialPreconditionerRuntime mTrivialPreconditionerRuntime;

private:
};

} // namespace lama

#endif // LAMA_TRIVIALPRECONDITIONER_HPP_
