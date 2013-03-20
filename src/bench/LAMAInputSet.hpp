/**
 * @file LAMAInputSet.hpp
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
 * @brief LAMAInputSet.hpp
 * @author Thomas Brandes, Jiri Kraus
 * @date 12.09.2011
 * $Id$
 */

#ifndef LAMA_LAMAINPUTSET_HPP_
#define LAMA_LAMAINPUTSET_HPP_

#include <memory>

#include <bench/config.hpp>

#include <framework/src/benchmark_framework.hpp>

#include <logging/logging.hpp>

#include <lama/matrix/CSRSparseMatrix.hpp>
#include <lama/matrix/DenseMatrix.hpp>

#include <lama/DenseVector.hpp>
#include <lama/NonCopyable.hpp>

/** Input set for LAMA is always distributed CSR sparse matrix with two vectors,
 all of type double.
 */

class LAMABENCH_DLL_IMPORTEXPORT LAMAInputSet: public bf::InputSet, private lama::NonCopyable
{
public:

    LAMAInputSet( const std::string& id );

    LAMAInputSet(
        const std::string& id,
        double alpha,
        double beta,
        std::auto_ptr<lama::DenseVector<double> > x,
        std::auto_ptr<lama::CSRSparseMatrix<double> > A );

    LAMAInputSet(
        const std::string& id,
        double alpha,
        double beta,
        std::auto_ptr<lama::DenseVector<double> > x,
        std::auto_ptr<lama::DenseVector<double> > y,
        std::auto_ptr<lama::CSRSparseMatrix<double> > A );

    LAMAInputSet(
        const std::string& id,
        double alpha,
        double beta,
        std::auto_ptr<lama::DenseVector<double> > x,
        std::auto_ptr<lama::DenseVector<double> > y,
        std::auto_ptr<lama::DenseMatrix<double> > A );

    virtual ~LAMAInputSet();

    double getAlpha() const;
    double getBeta() const;
    const lama::DenseVector<double>& getX() const;
    const lama::DenseVector<double>& getY() const;
    const lama::CSRSparseMatrix<double>& getA() const;
    const lama::CSRSparseMatrix<double>& getB() const;
    const lama::CSRSparseMatrix<double>& getC() const;

    const lama::DenseMatrix<double>& getDenseA() const;
    const lama::DenseMatrix<double>& getDenseB() const;
    const lama::DenseMatrix<double>& getDenseC() const;

private:

    LAMA_LOG_DECL_STATIC_LOGGER(logger);

    double mAlpha;
    double mBeta;

    lama::DenseVector<double>* mX;
    lama::DenseVector<double>* mY;

    mutable lama::CSRSparseMatrix<double>* mA;
    mutable lama::CSRSparseMatrix<double>* mB;
    mutable lama::CSRSparseMatrix<double>* mC;

    mutable lama::DenseMatrix<double>* mDenseA;
    mutable lama::DenseMatrix<double>* mDenseB;
    mutable lama::DenseMatrix<double>* mDenseC;
};

#endif // LAMA_LAMAINPUTSET_HPP_
