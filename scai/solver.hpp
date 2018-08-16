/**
 * @file solver.hpp
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
 * @brief General header file that includes all solver header files.
 * @author Lauretta Schubert
 * @date 15.01.2016
 */

#pragma once

#include <scai/solver/BiCGstab.hpp>
#include <scai/solver/CG.hpp>
#include <scai/solver/CGS.hpp>
#include <scai/solver/GMRES.hpp>
#include <scai/solver/InverseSolver.hpp>
#include <scai/solver/Jacobi.hpp>
#include <scai/solver/MINRES.hpp>
#include <scai/solver/QMR.hpp>
#include <scai/solver/Richardson.hpp>
#include <scai/solver/SimpleAMG.hpp>
#include <scai/solver/TFQMR.hpp>
#include <scai/solver/TrivialPreconditioner.hpp>

#include <scai/solver/logger/CommonLogger.hpp>

#include <scai/solver/criteria/IterationCount.hpp>
#include <scai/solver/criteria/ResidualThreshold.hpp>

