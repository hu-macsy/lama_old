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

#pragma once

#include <scai/common/unique_ptr.hpp>

#include <scai/benchmark.hpp>

#include <scai/logging.hpp>

#include <scai/lama/matrix/CSRSparseMatrix.hpp>
#include <scai/lama/matrix/DenseMatrix.hpp>

#include <scai/lama/DenseVector.hpp>
#include <scai/common/NonCopyable.hpp>

/** Input set for LAMA is always distributed CSR sparse matrix with two vectors,
 all of type double.
 */

class COMMON_DLL_IMPORTEXPORT LAMAInputSet: public bf::InputSet, private scai::common::NonCopyable
{
public:

    LAMAInputSet( const std::string& id );

    LAMAInputSet(
        const std::string& id,
        double alpha,
        double beta,
        scai::common::unique_ptr<scai::lama::DenseVector<double> > x,
        scai::common::unique_ptr<scai::lama::CSRSparseMatrix<double> > A );

    LAMAInputSet(
        const std::string& id,
        double alpha,
        double beta,
        scai::common::unique_ptr<scai::lama::DenseVector<double> > x,
        scai::common::unique_ptr<scai::lama::DenseVector<double> > y,
        scai::common::unique_ptr<scai::lama::CSRSparseMatrix<double> > A );

    LAMAInputSet(
        const std::string& id,
        double alpha,
        double beta,
        scai::common::unique_ptr<scai::lama::DenseVector<double> > x,
        scai::common::unique_ptr<scai::lama::DenseVector<double> > y,
        scai::common::unique_ptr<scai::lama::DenseMatrix<double> > A );

    virtual ~LAMAInputSet();

    double getAlpha() const;
    double getBeta() const;
    const scai::lama::DenseVector<double>& getX() const;
    const scai::lama::DenseVector<double>& getY() const;
    const scai::lama::CSRSparseMatrix<double>& getA() const;
    const scai::lama::CSRSparseMatrix<double>& getB() const;
    const scai::lama::CSRSparseMatrix<double>& getC() const;

    const scai::lama::DenseMatrix<double>& getDenseA() const;
    const scai::lama::DenseMatrix<double>& getDenseB() const;
    const scai::lama::DenseMatrix<double>& getDenseC() const;

private:

    SCAI_LOG_DECL_STATIC_LOGGER(logger);

    double mAlpha;
    double mBeta;

    scai::lama::DenseVector<double>* mX;
    scai::lama::DenseVector<double>* mY;

    mutable scai::lama::CSRSparseMatrix<double>* mA;
    mutable scai::lama::CSRSparseMatrix<double>* mB;
    mutable scai::lama::CSRSparseMatrix<double>* mC;

    mutable scai::lama::DenseMatrix<double>* mDenseA;
    mutable scai::lama::DenseMatrix<double>* mDenseB;
    mutable scai::lama::DenseMatrix<double>* mDenseC;
};
