/**
 * @file PMVBenchmark.cpp
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
 * @brief PMVBenchmark.cpp
 * @author brandes
 * @date 11.05.2011
 */

#include <scai/lama/benchmark/PMVBenchmark.hpp>

#include <scai/common/macros/loop.hpp>

// Matrix storage types needed for instantiation

#include <scai/lama/matrix/CSRSparseMatrix.hpp>
#include <scai/lama/matrix/ELLSparseMatrix.hpp>
#include <scai/lama/matrix/DIASparseMatrix.hpp>
#include <scai/lama/matrix/JDSSparseMatrix.hpp>
#include <scai/lama/matrix/COOSparseMatrix.hpp>
#include <scai/lama/matrix/DenseMatrix.hpp>

using namespace lama;

#define INSTANTIATE_PMV( ValueType )                                  \
    template class PMVBenchmark<CSRSparseMatrix<ValueType> >;        \
    template class PMVBenchmark<ELLSparseMatrix<ValueType> >;        \
    template class PMVBenchmark<JDSSparseMatrix<ValueType> >;        \
    template class PMVBenchmark<COOSparseMatrix<ValueType> >;        \
    template class PMVBenchmark<DIASparseMatrix<ValueType> >;        \
    template class PMVBenchmark<DenseMatrix<ValueType> >;

SCAI_COMMON_LOOP( INSTANTIATE_PMV , SCAI_NUMERIC_TYPES_HOST )

