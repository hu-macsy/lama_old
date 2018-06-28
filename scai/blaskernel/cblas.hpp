/**
 * @file cblas.hpp
 *
 * @license
 * Copyright (c) 2009-2018
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the SCAI framework LAMA.
 *
 * LAMA is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * LAMA is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
 * more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with LAMA. If not, see <http://www.gnu.org/licenses/>.
 * @endlicense
 *
 * @brief C Interface to BLAS routines
 * @author Thomas Brandes
 * @date 05.06.2014
 */

#pragma once

enum CBLAS_ORDER
{   CblasRowMajor = 101, CblasColMajor = 102};
enum CBLAS_TRANSPOSE
{   CblasNoTrans = 111, CblasTrans = 112, CblasConjTrans = 113};
enum CBLAS_UPLO
{   CblasUpper = 121, CblasLower = 122};
enum CBLAS_DIAG
{   CblasNonUnit = 131, CblasUnit = 132};
enum CBLAS_SIDE
{   CblasLeft = 141, CblasRight = 142};
