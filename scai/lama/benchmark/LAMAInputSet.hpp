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

#include <scai/benchmark.hpp>

#include <scai/logging.hpp>

#include <scai/lama/matrix/CSRSparseMatrix.hpp>
#include <scai/lama/matrix/DenseMatrix.hpp>

#include <scai/lama/DenseVector.hpp>
#include <scai/common/NonCopyable.hpp>

/** Input set for LAMA is always distributed CSR sparse matrix with two vectors,
 *   all of type double.
 */

namespace scai
{

namespace lama
{

/** This kind of input set consists of scalars alpha and beta, vectors X and Y, and matrix A */

class COMMON_DLL_IMPORTEXPORT LAMAInputSet: public benchmark::InputSet
{
public:

    virtual ~LAMAInputSet();

    /** Override implementation InputSet::writeAt */

    virtual void writeAt( std::ostream& stream ) const;

    double getAlpha() const;
    double getBeta() const;

    const DenseVector<double>& getX() const;
    const DenseVector<double>& getY() const;

    const CSRSparseMatrix<double>& getA() const;

protected:

    LAMAInputSet();

    SCAI_LOG_DECL_STATIC_LOGGER( logger );

    double mAlpha;
    double mBeta;

    std::unique_ptr<DenseVector<double> > mX;
    std::unique_ptr<DenseVector<double> > mY;

    std::unique_ptr<CSRSparseMatrix<double> > mA;
};

}

}
