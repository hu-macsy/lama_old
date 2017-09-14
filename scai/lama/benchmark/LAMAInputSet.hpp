/**
 * @file LAMAInputSet.hpp
 *
 * @license
 * Copyright (c) 2009-2017
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the SCAI framework LAMA.
 *
 * LAMA is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Affero General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * LAMA is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for
 * more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with LAMA. If not, see <http://www.gnu.org/licenses/>.
 *
 * Other Usage
 * Alternatively, this file may be used in accordance with the terms and
 * conditions contained in a signed written agreement between you and
 * Fraunhofer SCAI. Please contact our distributor via info[at]scapos.com.
 * @endlicense
 *
 * @brief LAMAInputSet.hpp
 * @author Thomas Brandes, Jiri Kraus
 * @date 12.09.2011
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

    SCAI_LOG_DECL_STATIC_LOGGER( logger );

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
