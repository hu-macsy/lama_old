/**
 * @file PardisoError.cpp
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
 * @brief Contains the implementation of the class PardisoError.
 * @author Lauretta Schubert
 * @date 18.07.2016
 */

namespace scai
{

namespace sparsekernel
{

const char* pardisoErrorString( int error )
{
    const char* str = "";

    switch ( error )
    {
        case 0:
            str = "Pardiso successful";
            break;
        case -1:
            str = "Pardiso: input inconsistent.";
            break;
        case -2:
            str = "Pardiso: not enough memory.";
            break;
        case -3:
            str = "Pardiso: reordering problem.";
            break;
        case -4:
            str = "Pardiso: zero pivot, numerical factorization or iterative refinement problem.";
            break;
        case -5:
            str = "Pardiso: unclassified (internal) error.";
            break;
        case -6:
            str = "Pardiso: reordering failed (matrix types 11 and 13 only)";
            break;
        case -7:
            str = "Pardiso: diagonal matrix is singular.";
            break;
        case -8:
            str = "Pardiso: 32-bit integer overflow problem.";
            break;
        case -9:
            str = "Pardiso: not enough memory for OOC.";
            break;
        case -10:
            str = "Pardiso: error opening OOC files.";
            break;
        case -11:
            str = "Pardiso: read/write error with OOC files.";
            break;
        case -12:
            str = "Pardiso: (pardiso_64 only) pardiso_64 called from 32-bit library.";
            break;

        default:
            str = "Unknown Pardiso error";
    }

    return str;

}

} /* end namespace sparsekernel */

} /* end namespace scai */
