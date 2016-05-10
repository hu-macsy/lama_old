/**
 * @file solver.hpp
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
 * @brief General header file that includes all solver header files.
 * @author Lauretta Schubert
 * @date 15.01.2016
 */

#pragma once

#include <scai/solver/BiCG.hpp>
#include <scai/solver/BiCGstab.hpp>
#include <scai/solver/CG.hpp>
#include <scai/solver/CGNE.hpp>
#include <scai/solver/CGNR.hpp>
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
